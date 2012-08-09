% function create_fine_grid(ncfile_coarse,ncfile_fine,Istr,Iend,Jstr,Jend,scale)
%
% called from create_child_grid.m
%
% user provides: - netcdf file name of coarse grid
%                - netcdf file name for new grid
%                - scale factor of increased resolution
%
% output is a netcdf roms grid with finer resolution
%
% User will need to modify the masking and bathymetry of new grid!!!!!
%
% Version jcwarner Aug 3, 2007
% Updated clea denamiel Feb 5, 2012
% Updated jcwarner 26Jun2012, make this a main driver.

%get the coarser grid
netcdf_load(ncfile_coarse)

%set some arrays for interpolation here, psi points
% _c for coarse, _f for fine
offset=ceil(4/scale);  % need 4 points to the left, 3 to the right
numx_c=((Iend+1)-(Istr-offset))+1;
numy_c=((Jend+1)-(Jstr-offset))+1;
numx_f=scale*(numx_c-1)+1;
numy_f=scale*(numy_c-1)+1;

xc_c=[1:numx_c];
xc_c=repmat(xc_c',1,numy_c);
yc_c=[1:numy_c];
yc_c=repmat(yc_c,numx_c,1);

xc_f=[1:1/scale:numx_c];
xc_f=repmat(xc_f',1,numy_f);
yc_f=[1:1/scale:numy_c];
yc_f=repmat(yc_f,numx_f,1);

%%%%%%%%%  Lon Lat SPACE  %%%%%%%%%%%
%establish a full grid of psi points, including an extra row and column
lon_full_grid=interp2(xc_c', yc_c', lon_psi(Istr-offset:Iend+1,Jstr-offset:Jend+1)', xc_f, yc_f);
lat_full_grid=interp2(xc_c', yc_c', lat_psi(Istr-offset:Iend+1,Jstr-offset:Jend+1)', xc_f, yc_f);

%set the psi points
fine.lon.psi=lon_full_grid(offset+2:end-(scale-3+1),offset+2:end-(scale-3+1)); 
fine.lat.psi=lat_full_grid(offset+2:end-(scale-3+1),offset+2:end-(scale-3+1)); 

%set the rho points
ztemp=interp2(lon_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lon.rho=ztemp(2:2:end-1,2:2:end-1);
ztemp=interp2(lat_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lat.rho=ztemp(2:2:end-1,2:2:end-1);

%set the u points 
ztemp=interp2(lon_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1); 
fine.lon.u=ztemp(3:2:end-2,2:2:end-1); 
ztemp=interp2(lat_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1); 
fine.lat.u=ztemp(3:2:end-2,2:2:end-1); 

%set the v points 
ztemp=interp2(lon_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1); 
fine.lon.v=ztemp(2:2:end-1,3:2:end-2); 
ztemp=interp2(lat_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1); 
fine.lat.v=ztemp(2:2:end-1,3:2:end-2); 

%%%%%%%%%  X Y SPACE  %%%%%%%%%%%
%establish a full grid of psi points, including an extra row and column
x_full_grid=interp2(xc_c', yc_c', x_psi(Istr-offset:Iend+1,Jstr-offset:Jend+1)', xc_f, yc_f);
y_full_grid=interp2(xc_c', yc_c', y_psi(Istr-offset:Iend+1,Jstr-offset:Jend+1)', xc_f, yc_f);

%set the psi points
fine.x.psi=x_full_grid(offset+2:end-(scale-3+1),offset+2:end-(scale-3+1)); 
fine.y.psi=y_full_grid(offset+2:end-(scale-3+1),offset+2:end-(scale-3+1)); 

%set the rho points
ztemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.x.rho=ztemp(2:2:end-1,2:2:end-1);
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.y.rho=ztemp(2:2:end-1,2:2:end-1);

%set the u points 
ztemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1); 
fine.x.u=ztemp(3:2:end-2,2:2:end-1); 
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1); 
fine.y.u=ztemp(3:2:end-2,2:2:end-1); 

%set the v points 
ztemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1); 
fine.x.v=ztemp(2:2:end-1,3:2:end-2); 
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1); 
fine.y.v=ztemp(2:2:end-1,3:2:end-2); 

%%%%%%%  now need grid metrics %%%%%%%%%%%
fine.pm=zeros(size(fine.lon.rho));
fine.pn=zeros(size(fine.lat.rho));
for j=Jstr-1:Jend+1
  for i=Istr-1:Iend+1
    X=[lon_psi(i-1,j-1) lon_psi(i,j-1) lon_psi(i,j) lon_psi(i-1,j)];
    Y=[lat_psi(i-1,j-1) lat_psi(i,j-1) lat_psi(i,j) lat_psi(i-1,j)];
    %find child cells within parent cell
    ind=find(inpolygon(fine.lon.rho(:),fine.lat.rho(:),X,Y));
    if ~isempty(ind),
      fine.pm(ind)=pm(i,j)*scale;
      fine.pn(ind)=pn(i,j)*scale;
    end
  end
end

fine.dmde = zeros(size(fine.pm));
fine.dndx = zeros(size(fine.pn));
fine.dndx(2:end-1, :) = 0.5*(1./fine.pn(3:end, :) - 1./fine.pn(1:end-2, :));
fine.dmde(:, 2:end-1) = 0.5*(1./fine.pm(:, 3:end) - 1./fine.pm(:, 1:end-2));
fine.dmde(:,1)=fine.dmde(:,2);
fine.dmde(:,end)=fine.dmde(:,end-1);
fine.dndx(1,:)=fine.dndx(2,:);
fine.dndx(end,:)=fine.dndx(end-1,:);

% Grid-cell orientation, degrees counter-clockwise from 
fine.angle=griddata(lon_rho,lat_rho,angle,fine.lon.rho,fine.lat.rho);
fine.f=griddata(lon_rho,lat_rho,f,fine.lon.rho,fine.lat.rho);

%%%%%% GRID h and masking  %%%%%%%%%
[XI,YI]=meshgrid([min(x_psi(:)):min(1./pm(:))/2:max(x_psi(:))], ...
                 [min(y_psi(:)):min(1./pn(:))/2:max(y_psi(:))]);
ZI=griddata(x_rho,y_rho,h,XI,YI);
%ZI=interp2(x_rho,y_rho,h,XI,YI);
fine.h=interp2(XI,YI,ZI,fine.x.rho,fine.y.rho);

ZI=griddata(x_rho,y_rho,mask_rho,XI,YI);
%ZI=interp2(x_rho,y_rho,mask_rho,XI,YI);
fine.mask.rho=interp2(XI,YI,ZI,fine.x.rho,fine.y.rho,'nearest');
%force masking to = 0 if it is less than 1.
fine.mask.rho(fine.mask.rho>0.95)=1;
fine.mask.rho(fine.mask.rho<0.96)=0;

%%%%%%%%%%%%%%%%%%% finished creating grid %%%%%%%%%%%%%%%%

%call routine to create netcdf file
[LP, MP] = size(fine.h);
create_roms_netcdf_grid_file(ncfile_fine,LP,MP)

%now open that file and prep for writing.
nc=netcdf.open(ncfile_fine,'NC_WRITE');

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
if (exist('JPRJ','var'))
  projection=JPRJ;
else
  projection='';
end
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
netcdf.putVar(nc,v4,spherical)

xl = max(fine.x.rho(:)) - min(fine.x.rho(:));
el = max(fine.y.rho(:)) - min(fine.y.rho(:));
netcdf.putVar(nc,v1,xl)
netcdf.putVar(nc,v2,el)

% Depths at RHO points.
bathymetry = fine.h;
if ~isempty(bathymetry)
  netcdf.putVar(nc,v5,min(min(bathymetry)));
  netcdf.putVar(nc,v6,max(max(bathymetry)));
  netcdf.putVar(nc,v7,[0 0 0],[LP MP 1],bathymetry);
  netcdf.putVar(nc,v8,bathymetry);
end

% Coriolis
netcdf.putVar(nc,v9,fine.f)

% Locations.
netcdf.putVar(nc,v14,fine.x.rho);
netcdf.putVar(nc,v15,fine.y.rho);

netcdf.putVar(nc,v16,fine.x.psi);
netcdf.putVar(nc,v17,fine.y.psi);

netcdf.putVar(nc,v18,fine.x.u);
netcdf.putVar(nc,v19,fine.y.u);

netcdf.putVar(nc,v20,fine.x.v);
netcdf.putVar(nc,v21,fine.y.v);

netcdf.putVar(nc,v22,fine.lat.rho);
netcdf.putVar(nc,v23,fine.lon.rho);

netcdf.putVar(nc,v24,fine.lat.psi);
netcdf.putVar(nc,v25,fine.lon.psi);

netcdf.putVar(nc,v26,fine.lat.u);
netcdf.putVar(nc,v27,fine.lon.u);

netcdf.putVar(nc,v28,fine.lat.v);
netcdf.putVar(nc,v29,fine.lon.v);

netcdf.putVar(nc,v10,fine.pm);
netcdf.putVar(nc,v11,fine.pn);

netcdf.putVar(nc,v12,fine.dndx);
netcdf.putVar(nc,v13,fine.dmde);

% Masking.
mask = fine.mask.rho;

water = double(mask);
netcdf.putVar(nc,v30,water);

u_mask = water(1:end-1,:) & water(2:end,:);
netcdf.putVar(nc,v31,double(u_mask));

v_mask= water(:,1:end-1) & water(:,2:end);
netcdf.putVar(nc,v32,double(v_mask));

psi_mask= water(1:end-1,1:end-1) & water(1:end-1,2:end) & water(2:end,1:end-1) & water(2:end,2:end);
netcdf.putVar(nc,v33,double(psi_mask));

% Angle.
netcdf.putVar(nc,v34,fine.angle);  % Degrees.

netcdf.close(nc)


%to create swan grid then use:
  disp([' '])
  if (spherical)
    roms2swan(fine.lon.rho,fine.lat.rho,fine.h,fine.mask.rho);
  else
    roms2swan(fine.x.rho,fine.y.rho,fine.h,fine.mask.rho);
  end
  disp(['Created swan grid + bathy files: swan_bathy.bot and swan_coord.grd'])
  disp(['You should probably rename these so they dont get overwritten. '])

%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

% Plot new grid points inside old grid
for mm=1:2
  if (mm==1); netcdf_load(ncfile_coarse); end
  if (mm==2); netcdf_load(ncfile_fine); end
  if (spherical)
    xr=lon_rho;  yr=lat_rho;
    xp=lon_psi;  yp=lat_psi;
    xu=lon_u;    yu=lat_u;
    xv=lon_v;    yv=lat_v;
  else     
    xr=x_rho;  yr=y_rho;
    xp=x_psi;  yp=y_psi;
    xu=x_u;    yu=y_u;
    xv=x_v;    yv=y_v;
  end
  maskr=mask_rho;
  maskr(mask_rho==0)=nan;
  
% Grid lines
  figure(1)
  if (mm==1) 
    plot(xp,yp,'c')
    hold on
    plot(xp',yp','c')
    x_rhoc=xr; y_rhoc=yr;  %save these to plot later
  else
    plot(xp,yp,'k')
    plot(xp',yp','k')
    plot(xr,yr,'ro')
    plot(xu,yu,'kx')
    plot(xv,yv,'ms')
    plot(x_rhoc,y_rhoc,'bo') %replot coarse rhos points on top
  end
  title('grid lines and rho, u, v, points')
% Bathymetry
  figure(2)
  pcolorjw(xr,yr,h.*maskr)
  hold on
  shading flat
  colorbar
  if (mm==2)
    title('bathymetry (m)')
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end
  
% PM
  figure(3)
  pcolorjw(xr,yr,pm.*maskr)
  shading flat
  hold on
  colorbar
  if (mm==2)
    title('pm');
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end

% PN
  figure(4)
  pcolorjw(xr,yr,pn.*maskr)
  shading flat
  hold on
  colorbar
  if (mm==2)
    title('pn');
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end

% dndx
  figure(5)
  pcolorjw(xr,yr,dndx.*maskr)
  shading flat
  hold on
  colorbar
  if (mm==2)
    title('dndx');
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end

% dmde
  figure(6)
  pcolorjw(xr,yr,dmde.*maskr)
  shading flat
  colorbar
  hold on
  if (mm==2)
    title('dmde');
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end

% angle
  figure(7)
  pcolorjw(xr,yr,angle.*maskr)
  shading flat
  hold on
  colorbar
  if (mm==2);
    title('angle');
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end
end

