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

%
% Call routine to create netcdf grid file and
% fill it with vars and attributes.
%
[Lp, Mp] = size(s.rho.depth);
L = Lp-1;   Lm = L-1;
M = Mp-1;   Mm = M-1;
create_roms_netcdf_grid_file(theROMSFile,Lp,Mp)

%now open that file and prep for writing.
nc=netcdf.open(theROMSFile,'NC_WRITE');

v1 = netcdf.inqVarID(nc,'xl');
v2 = netcdf.inqVarID(nc,'el');
v3 = netcdf.inqVarID(nc,'JPRJ');
v4 = netcdf.inqVarID(nc,'spherical');
v5 = netcdf.inqVarID(nc,'depthmin');
v6 = netcdf.inqVarID(nc,'depthmax');
v7 = netcdf.inqVarID(nc,'hraw');
v8 = netcdf.inqVarID(nc,'h');
v9 = netcdf.inqVarID(nc,'f');
v10 = netcdf.inqVarID(nc,'pm');
v11 = netcdf.inqVarID(nc,'pn');
v12 = netcdf.inqVarID(nc,'dndx');
v13 = netcdf.inqVarID(nc,'dmde');
v14 = netcdf.inqVarID(nc,'x_rho');
v15 = netcdf.inqVarID(nc,'y_rho');
v16 = netcdf.inqVarID(nc,'x_psi');
v17 = netcdf.inqVarID(nc,'y_psi');
v18 = netcdf.inqVarID(nc,'x_u');
v19 = netcdf.inqVarID(nc,'y_u');
v20 = netcdf.inqVarID(nc,'x_v');
v21 = netcdf.inqVarID(nc,'y_v');
v22 = netcdf.inqVarID(nc,'lat_rho');
v23 = netcdf.inqVarID(nc,'lon_rho');
v24 = netcdf.inqVarID(nc,'lat_psi');
v25 = netcdf.inqVarID(nc,'lon_psi');
v26 = netcdf.inqVarID(nc,'lat_u');
v27 = netcdf.inqVarID(nc,'lon_u');
v28 = netcdf.inqVarID(nc,'lat_v');
v29 = netcdf.inqVarID(nc,'lon_v');
v30 = netcdf.inqVarID(nc,'mask_rho');
v31 = netcdf.inqVarID(nc,'mask_u');
v32 = netcdf.inqVarID(nc,'mask_v');
v33 = netcdf.inqVarID(nc,'mask_psi');
v34 = netcdf.inqVarID(nc,'angle');

% Fill the variables with data.

disp(' ## Filling Variables...')

projection = s.projection;
m_proj(s.projection);
switch lower(projection)
%m_proj(s.projection.name);
%switch lower(projection.name)

case 'mercator'
	theProjection = 'ME';
case 'stereographic'
	theProjection = 'ST';
case 'lambert conformal conic'
	theProjection = 'LC';
    m_proj('Lambert Conformal Conic','clo',mean(s.rho.lon(:)))
case 'utm'
	theProjection = 'TM';
    m_proj('UTM','zone',s.projection.zone,'longitudes',[s.projection.lon1 s.projection.lon2]);
otherwise
	theProjection = '??';
end
netcdf.putVar(nc,v3,theProjection)
netcdf.putVar(nc,v4,s.spherical)

% Determine coordinates.
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
  netcdf.putVar(nc,v7,[0 0 0],[Lp Mp 1],bathymetry);
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
if (1)  
  if ((s.spherical=='T') || (s.spherical=='t'))
    lat_u=geogrid_lat(i_u, j_u);
    lon_u=geogrid_lon(i_u, j_u);
    [dx,ang]=sw_dist(lat_u(:),lon_u(:),'km');
  
    dx=[dx(:); dx(end)]*1000;  % km==> m
    dx=reshape(dx,Lp-1,Mp);
    dx=[dx(1,:); dx(1:end-1,:); dx(end-1,:)];
  
    ang=[ang(:); ang(end)];
    ang=reshape(ang,Lp-1,Mp);
    ang=[ang(1,:); ang(1:end-1,:); ang(end-1,:)];
    ang=ang*pi/180;
  
    lat_v=geogrid_lat(i_v, j_v).';
    lon_v=geogrid_lon(i_v, j_v).';
    dy=sw_dist(lat_v(:),lon_v(:),'km');
  
    dy=[dy(:); dy(end)]*1000;   % km ==> m 
    dy=reshape(dy,Mp-1,Lp);
    dy=[dy(1,:); dy(1:end-1,:); dy(end-1,:)];
    dy=dy.';
  else
    x_r=grid_x(i_rho, j_rho);
    y_r=grid_y(i_rho, j_rho);
  %
    x_u=grid_x(i_u, j_u);
    y_u=grid_y(i_u, j_u);
    %dx=diff(x_u);
    dx=sqrt((x_u(2:end,:)-x_u(1:end-1,:)).^2+(y_u(2:end,:)-y_u(1:end-1,:)).^2);
    dxl=sqrt((x_r(1,:)-x_u(1,:)).^2+(y_r(1,:)-y_u(1,:)).^2);
    dxr=sqrt((x_r(Lp,:)-x_u(L,:)).^2+(y_r(Lp,:)-y_u(L,:)).^2);
    dx=[2*dxl; dx; 2*dxr];
  
    x_v=grid_x(i_v, j_v);
    y_v=grid_y(i_v, j_v);
  % dy=diff(y_v);
    dy=sqrt((x_v(:,2:end)-x_v(:,1:end-1)).^2+(y_v(:,2:end)-y_v(:,1:end-1)).^2);
    dyl=sqrt((x_r(:,1)-x_v(:,1)).^2+(y_r(:,1)-y_v(:,1)).^2);
    dyr=sqrt((x_r(:,Mp)-x_v(:,M)).^2+(y_r(:,Mp)-y_v(:,M)).^2);
    dy=[2*dyl dy 2*dyr];
  
    y_v=grid_y(i_psi, j_psi);
    y_v=diff(y_v); 
    x_v=dx(2:end-1,2:end);
    if (isfield(s,'roms_angle'))
      ang=s.roms_angle;
    else
%     ang=angle(x_v+y_v*sqrt(-1));
%     ang=[ang(1,:); ang; ang(end,:)];
%     ang=[ang ang(:,end)];
      [~,dxdxi]=gradient(x);
      [~,dydxi]=gradient(y);
      ang=atan2(dydxi,dxdxi);
      ang(1,:)=ang(2,:)-(ang(3,:)-ang(2,:));
      ang(end,:)=ang(end-1,:)+(ang(end-1,:)-ang(end-2,:));
    end
  end
end

if (0)
  deg2rad = pi/180.0;  

  % Compute grid spacing.
  
  if ((s.spherical=='T') || (s.spherical=='t'))
      Xr=geogrid_lon(i_rho, j_rho);
      Yr=geogrid_lat(i_rho, j_rho);
      Xu=geogrid_lon(i_u, j_u);
      Yu=geogrid_lat(i_u, j_u);
      Xv=geogrid_lon(i_v, j_v); 
      Yv=geogrid_lat(i_v, j_v);
      dx = zeros(size(Xr));
      dy = zeros(size(Xr));
%
      x=(Xr-mean(Xr(:)))*deg2rad.*cosd(Yr);
      y=(Yr-mean(Yr(:)))*deg2rad;
      dx(2:L,1:Mp) = gcircle(Xu(1:Lm,1:Mp), Yu(1:Lm,1:Mp),                ...
          Xu(2:L ,1:Mp), Yu(2:L ,1:Mp));
      dx(1  ,1:Mp) = gcircle(Xr(1   ,1:Mp), Yr(1   ,1:Mp),                ...
          Xu(1   ,1:Mp), Yu(1   ,1:Mp)).*2.0;
      dx(Lp ,1:Mp) = gcircle(Xr(Lp  ,1:Mp), Yr(Lp  ,1:Mp),                ...
          Xu(L   ,1:Mp), Yu(L   ,1:Mp)).*2.0;
  
      dy(1:Lp,2:M) = gcircle(Xv(1:Lp,1:Mm), Yv(1:Lp,1:Mm),                ...
          Xv(1:Lp,2:M ), Yv(1:Lp,2:M ));
      dy(1:Lp,1  ) = gcircle(Xr(1:Lp,1   ), Yr(1:Lp,1   ),                ...
          Xv(1:Lp,1   ), Yv(1:Lp,1   )).*2.0;
      dy(1:Lp,Mp ) = gcircle(Xr(1:Lp,Mp  ), Yr(1:Lp,Mp  ),                ...
          Xv(1:Lp,M   ), Yv(1:Lp,M   )).*2.0;
  
      dx = dx .* 1000;       % great circle function computes
      dy = dy .* 1000;       % distances in kilometers
  else                                      % Cartesian distance
      Xr=grid_x(i_rho, j_rho);
      Yr=grid_y(i_rho, j_rho);
      Xu=grid_x(i_u, j_u);
      Yu=grid_y(i_u, j_u);
      Xv=grid_x(i_v, j_v);
      Yv=grid_y(i_v, j_v);
      dx = zeros(size(Xr));
      dy = zeros(size(Xr));
%
      x=Xr;
      y=Yr;
      dx(2:L,1:Mp) = sqrt((Xu(2:L ,1:Mp) - Xu(1:Lm,1:Mp)).^2 +            ...
          (Yu(2:L ,1:Mp) - Yu(1:Lm,1:Mp)).^2);
      dx(1  ,1:Mp) = sqrt((Xu(1   ,1:Mp) - Xr(1   ,1:Mp)).^2 +            ...
          (Yu(1   ,1:Mp) - Yr(1   ,1:Mp)).^2).*2.0;
      dx(Lp ,1:Mp) = sqrt((Xr(Lp  ,1:Mp) - Xu(L   ,1:Mp)).^2 +            ...
          (Yr(Lp  ,1:Mp) - Yu(L   ,1:Mp)).^2).*2.0;
  
      dy(1:Lp,2:M) = sqrt((Xv(1:Lp,2:M ) - Xv(1:Lp,1:Mm)).^2 +            ...
          (Yv(1:Lp,2:M ) - Yv(1:Lp,1:Mm)).^2);
      dy(1:Lp,1  ) = sqrt((Xv(1:Lp,1   ) - Xr(1:Lp,1   )).^2 +            ...
          (Yv(1:Lp,1   ) - Yr(1:Lp,1   )).^2).*2.0;
      dy(1:Lp,Mp ) = sqrt((Xr(1:Lp,Mp  ) - Xv(1:Lp,M   )).^2 +            ...
          (Yr(1:Lp,Mp  ) - Yv(1:Lp,M   )).^2).*2.0;
  end
%
  [~,dxdxi]=gradient(x);
  [~,dydxi]=gradient(y);
  ang=atan2(dydxi,dxdxi);
  ang(1,:)=ang(2,:)-(ang(3,:)-ang(2,:));
  ang(end,:)=ang(end-1,:)+(ang(end-1,:)-ang(end-2,:));
end

pm=1./dx;
pn=1./dy;
netcdf.putVar(nc,v10,pm);
netcdf.putVar(nc,v11,pn);

dmde = zeros(Lp, Mp);
dndx = zeros(Lp, Mp);
dmde(:,2:end-1) = 0.5*(dx(:,3:end) - dx(:,1:end-2));
dmde(:,1)=dmde(:,2);
dmde(:,end)=dmde(:,end-1);
dndx(2:end-1,:) = 0.5*(dy(3:end,:) - dy(1:end-2,:));
dndx(1,:)=dndx(2,:);
dndx(end,:)=dndx(end-1,:);
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
