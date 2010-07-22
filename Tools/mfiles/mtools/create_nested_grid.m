% create_nested_grid.m
%
% user provides: - netcdf file name of coarse grid
%                - netcdf file name for new grid
%                - scale factor of increased resolution
%
% output is a netcdf roms grid with finer resolution
%
% User will need to modify the masking and bathymetry of new grid!!!!!
%
% jcwarner Aug 3, 2007
%

%%%%%%%%%%%%% START OF USER SECTION %%%%%%%%%%%%%%%%

%1) ENTER NAME OF EXISITNG COARSE GRID FILE
ncfile_coarse='D:\data\models\COAWST\Projects\Inlet_test\Refined\inlet_test_grid.nc';

%2) ENTER NAME OF NEW FINE GRID FILE
ncfile_fine='inlet_test_grid_ref5.nc';

%3) ENTER THE START AND END INDICES OF THE PSI POINTS OF THE COARSE GRID
%   THAT IDENTIFY THE OUTER BOUNDS OF THE FINE GRID
Istr=24;
Iend=54;
Jstr=40;
Jend=56;

%4) ENTER SCALE FACTOR FOR INCREASED RESOLUTION (use 3 or 5)
scale=5;

%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

%set some random variables here
msize1=5;      %marker size for plots

%get the coarser grid
ncload(ncfile_coarse)

%size of coarse grid section to have increased resolution
Ilen_c=Iend-Istr;
Jlen_c=Jend-Jstr;

%xi direction (left-right)
LP = Ilen_c*scale+2;   %rho dimension
L = LP-1;              %psi dimension

%eta direction (up-down)
MP = Jlen_c*scale+2;   % rho dimension
M = MP-1;              % psi dimension

%sizes for fine grid
xi_psi = L;
xi_rho = LP;
xi_u = L;
xi_v = LP;

eta_psi = M;
eta_rho = MP;
eta_u = MP;
eta_v = M;

%set some arrays for interpolation here
% _c for coarse, _f for fine
xc_c=[1:scale:xi_psi];
xc_c=repmat(xc_c,Jlen_c+1,1);
yc_c=[1:scale:eta_psi];
yc_c=repmat(yc_c',1,Ilen_c+1);

xc_f=[1:xi_psi];
xc_f=repmat(xc_f,eta_psi,1);
yc_f=[1:eta_psi];
yc_f=repmat(yc_f',1,xi_psi);


%%%%%%%%%  Lon Lat SPACE  %%%%%%%%%%%
%interpolate the psi points
fine.lon.psi=interp2(xc_c, yc_c, lon_psi(Jstr:Jend,Istr:Iend), xc_f, yc_f);
fine.lat.psi=interp2(xc_c, yc_c, lat_psi(Jstr:Jend,Istr:Iend), xc_f, yc_f);

%now make up a fake set of psi points that includes one extra bounding
%row and column on all sides.
  x_full_grid=ones(size(fine.lon.psi)+2);
  x_full_grid(2:end-1,2:end-1)=fine.lon.psi;
  x_full_grid(2:end-1,1)=fine.lon.psi(:,1)-(fine.lon.psi(:,2)-fine.lon.psi(:,1));
  x_full_grid(2:end-1,end)=fine.lon.psi(:,end)+(fine.lon.psi(:,end)-fine.lon.psi(:,end-1));
  x_full_grid(1,2:end-1)=fine.lon.psi(1,:)-(fine.lon.psi(2,:)-fine.lon.psi(1,:));
  x_full_grid(end,2:end-1)=fine.lon.psi(end,:)+(fine.lon.psi(end,:)-fine.lon.psi(end-1,:));
  x_full_grid(1,1)=x_full_grid(1,2)-(x_full_grid(1,3)-x_full_grid(1,2));
  x_full_grid(1,end)=x_full_grid(1,end-1)+(x_full_grid(1,end-1)-x_full_grid(1,end-2));
  x_full_grid(end,1)=x_full_grid(end,2)-(x_full_grid(end,3)-x_full_grid(end,2));
  x_full_grid(end,end)=x_full_grid(end,end-1)+(x_full_grid(end,end-1)-x_full_grid(end,end-2));
%
  y_full_grid=ones(size(fine.lat.psi)+2);
  y_full_grid(2:end-1,2:end-1)=fine.lat.psi;
  y_full_grid(2:end-1,1)=fine.lat.psi(:,1)-(fine.lat.psi(:,2)-fine.lat.psi(:,1));
  y_full_grid(2:end-1,end)=fine.lat.psi(:,end)+(fine.lat.psi(:,end)-fine.lat.psi(:,end-1));
  y_full_grid(1,2:end-1)=fine.lat.psi(1,:)-(fine.lat.psi(2,:)-fine.lat.psi(1,:));
  y_full_grid(end,2:end-1)=fine.lat.psi(end,:)+(fine.lat.psi(end,:)-fine.lat.psi(end-1,:));
  y_full_grid(1,1)=y_full_grid(1,2)-(y_full_grid(1,3)-y_full_grid(1,2));
  y_full_grid(1,end)=y_full_grid(1,end-1)+(y_full_grid(1,end-1)-y_full_grid(1,end-2));
  y_full_grid(end,1)=y_full_grid(end,2)-(y_full_grid(end,3)-y_full_grid(end,2));
  y_full_grid(end,end)=y_full_grid(end,end-1)+(y_full_grid(end,end-1)-y_full_grid(end,end-2));

%set the rho points
ztemp=interp2(x_full_grid,1);
fine.lon.rho=ztemp(2:2:end-1,2:2:end-1);
ztemp=interp2(y_full_grid,1);
fine.lat.rho=ztemp(2:2:end-1,2:2:end-1);

%set the u points
ztemp=interp2(x_full_grid,1);
fine.lon.u=ztemp(2:2:end-1,3:2:end-2);
ztemp=interp2(y_full_grid,1);
fine.lat.u=ztemp(2:2:end-1,3:2:end-2);

%set the v points
ztemp=interp2(x_full_grid,1);
fine.lon.v=ztemp(3:2:end-2,2:2:end-1);
ztemp=interp2(y_full_grid,1);
fine.lat.v=ztemp(3:2:end-2,2:2:end-1);


%%%%%%%%%  X Y SPACE  %%%%%%%%%%%
%interpolate the psi points
fine.x.psi=interp2(xc_c, yc_c, x_psi(Jstr:Jend,Istr:Iend), xc_f, yc_f);
fine.y.psi=interp2(xc_c, yc_c, y_psi(Jstr:Jend,Istr:Iend), xc_f, yc_f);

%now make up a fake set of psi points that includes one extra bounding
%row and column on all sides.
  x_full_grid=ones(size(fine.x.psi)+2);
  x_full_grid(2:end-1,2:end-1)=fine.x.psi;
  x_full_grid(2:end-1,1)=fine.x.psi(:,1)-(fine.x.psi(:,2)-fine.x.psi(:,1));
  x_full_grid(2:end-1,end)=fine.x.psi(:,end)+(fine.x.psi(:,end)-fine.x.psi(:,end-1));
  x_full_grid(1,2:end-1)=fine.x.psi(1,:)-(fine.x.psi(2,:)-fine.x.psi(1,:));
  x_full_grid(end,2:end-1)=fine.x.psi(end,:)+(fine.x.psi(end,:)-fine.x.psi(end-1,:));
  x_full_grid(1,1)=x_full_grid(1,2)-(x_full_grid(1,3)-x_full_grid(1,2));
  x_full_grid(1,end)=x_full_grid(1,end-1)+(x_full_grid(1,end-1)-x_full_grid(1,end-2));
  x_full_grid(end,1)=x_full_grid(end,2)-(x_full_grid(end,3)-x_full_grid(end,2));
  x_full_grid(end,end)=x_full_grid(end,end-1)+(x_full_grid(end,end-1)-x_full_grid(end,end-2));
%
  y_full_grid=ones(size(fine.y.psi)+2);
  y_full_grid(2:end-1,2:end-1)=fine.y.psi;
  y_full_grid(2:end-1,1)=fine.y.psi(:,1)-(fine.y.psi(:,2)-fine.y.psi(:,1));
  y_full_grid(2:end-1,end)=fine.y.psi(:,end)+(fine.y.psi(:,end)-fine.y.psi(:,end-1));
  y_full_grid(1,2:end-1)=fine.y.psi(1,:)-(fine.y.psi(2,:)-fine.y.psi(1,:));
  y_full_grid(end,2:end-1)=fine.y.psi(end,:)+(fine.y.psi(end,:)-fine.y.psi(end-1,:));
  y_full_grid(1,1)=y_full_grid(1,2)-(y_full_grid(1,3)-y_full_grid(1,2));
  y_full_grid(1,end)=y_full_grid(1,end-1)+(y_full_grid(1,end-1)-y_full_grid(1,end-2));
  y_full_grid(end,1)=y_full_grid(end,2)-(y_full_grid(end,3)-y_full_grid(end,2));
  y_full_grid(end,end)=y_full_grid(end,end-1)+(y_full_grid(end,end-1)-y_full_grid(end,end-2));

%set the rho points
ztemp=interp2(x_full_grid,1);
fine.x.rho=ztemp(2:2:end-1,2:2:end-1);
ztemp=interp2(y_full_grid,1);
fine.y.rho=ztemp(2:2:end-1,2:2:end-1);

%set the u points
ztemp=interp2(x_full_grid,1);
fine.x.u=ztemp(2:2:end-1,3:2:end-2);
ztemp=interp2(y_full_grid,1);
fine.y.u=ztemp(2:2:end-1,3:2:end-2);

%set the v points
ztemp=interp2(x_full_grid,1);
fine.x.v=ztemp(3:2:end-2,2:2:end-1);
ztemp=interp2(y_full_grid,1);
fine.y.v=ztemp(3:2:end-2,2:2:end-1);

%%%%%% GRID Metrics  %%%%%%%%%

[XI,YI]=meshgrid([min(x_psi(:)):min(1./pm(:))/2:max(x_psi(:))], ...
                 [min(y_psi(:)):min(1./pn(:))/2:max(y_psi(:))]);
%ZI=griddata(x_rho,y_rho,h,XI,YI);
ZI=interp2(x_rho,y_rho,h,XI,YI);
fine.h=interp2(XI,YI,ZI,fine.x.rho,fine.y.rho);
ZI=griddata(x_rho,y_rho,f,XI,YI);
%ZI=interp2(x_rho,y_rho,f,XI,YI);
fine.f=interp2(XI,YI,ZI,fine.x.rho,fine.y.rho);
ZI=griddata(x_rho,y_rho,mask_rho,XI,YI);
%ZI=interp2(x_rho,y_rho,mask_rho,XI,YI);
fine.mask.rho=interp2(XI,YI,ZI,fine.x.rho,fine.y.rho);

%force masking to = 0 if it is less than 1.
[js,is]=size(fine.mask.rho);
for i=1:is
  for j=1:js
    if (fine.mask.rho(j,i) < 1)
      fine.mask.rho(j,i)=0;
    end
  end
end
for i = 2:LP
  for j = 1:MP
    fine.mask.u(j, i-1) = fine.mask.rho(j, i) * fine.mask.rho(j, i-1);
  end
end
for i = 1:LP
  for j = 2:MP
    fine.mask.v(j-1, i) = fine.mask.rho(j, i) * fine.mask.rho(j-1, i);
  end
end
for i = 2:LP
  for j = 2:MP
    fine.mask.psi(j-1, i-1) = fine.mask.rho(j, i) * fine.mask.rho(j, i-1) * ...
                              fine.mask.rho(j-1, i) *fine.mask.rho(j-1, i-1);
  end
end

% compute dx and dy at u and v points for future computation of
% pm, pn, dmde, dndx
for i=1:LP
  for j=1:MP+1
    dx_v(j,i)=x_full_grid(j,i+1)-x_full_grid(j,i);
  end
end
for i=1:LP+1
  for j=1:MP
    dy_u(j,i)=y_full_grid(j+1,i)-y_full_grid(j,i);
  end
end
for i=1:LP
  for j=1:MP+1
    dy_v(j,i)=y_full_grid(j,i+1)-y_full_grid(j,i);
  end
end

%compute pm and pn
for i=1:LP
  for j=1:MP
    fine.pm(j,i)=1/(0.5*(dx_v(j,i)+dx_v(j+1,i  )));
    fine.pn(j,i)=1/(0.5*(dy_u(j,i)+dy_u(j  ,i+1)));
  end
end

%compute dmde and dndx
for i=1:LP
  for j=1:MP
    fine.dmde(j,i)=dx_v(j+1,i)-dx_v(j,i);
    fine.dndx(j,i)=dy_u(j,i+1)-dy_u(j,i);
  end
end

% compute angle
for i=1:LP
  for j=1:MP+1
    angle_v(j,i)=atan2(dy_v(j,i),dx_v(j,i));
  end
end
for i=1:LP
  for j=1:MP
    fine.angle(j,i)=0.5*(angle_v(j+1,i)+angle_v(j,i));
  end
end

%%%%%%%%%%%%%%%%%%% finished creating grid %%%%%%%%%%%%%%%%

%call routine to create netcdf file
roms_grid_create_netcdf(ncfile_fine,LP,MP)

%now fill that netcdf file
eval(['nc=netcdf(''',ncfile_fine,''',''w'');'])

xl = max(fine.x.psi(:)) - min(fine.x.psi(:));
el = max(fine.y.psi(:)) - min(fine.y.psi(:));
nc{'xl'}(:) = xl;
nc{'el'}(:) = el;

nc{'f'}(:) = fine.f;

% Locations.

nc{'x_rho'}(:) = fine.x.rho;
nc{'y_rho'}(:) = fine.y.rho;

nc{'x_psi'}(:) = fine.x.psi;
nc{'y_psi'}(:) = fine.y.psi;

nc{'x_u'}(:) = fine.x.u;
nc{'y_u'}(:) = fine.y.u;

nc{'x_v'}(:) = fine.x.v;
nc{'y_v'}(:) = fine.y.v;

nc{'lon_rho'}(:) = fine.lon.rho;
nc{'lat_rho'}(:) = fine.lat.rho;

nc{'lon_psi'}(:) = fine.lon.psi;
nc{'lat_psi'}(:) = fine.lat.psi;

nc{'lon_u'}(:) = fine.lon.u;
nc{'lat_u'}(:) = fine.lat.u;

nc{'lon_v'}(:) = fine.lon.v;
nc{'lat_v'}(:) = fine.lat.v;

% Metric factors.
nc{'dmde'}(:) = fine.dmde;
nc{'dndx'}(:) = fine.dndx;

nc{'pm'}(:) = fine.pm;
nc{'pn'}(:) = fine.pn;

% Angle.
nc{'angle'}(:) = fine.angle;

% Depths at RHO points.
nc{'h'}(:) = fine.h;
nc{'depthmin'}(:) = min(fine.h(:));
nc{'depthmax'}(:) = max(fine.h(:));

% Maksing. water = 1, land= 0.
nc{'mask_rho'}(:) = fine.mask.rho;
nc{'mask_u'}(:) = fine.mask.u;
nc{'mask_v'}(:) = fine.mask.v;
nc{'mask_psi'}(:) = fine.mask.psi;

% Close the ROMS File.

if ~isempty(close(nc))
	disp(' ## Unable too close the ROMS output file.')
end

ncclose('nc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot new grid inside old grid

figure
plot(x_psi,y_psi,'c')
hold on
plot(x_psi',y_psi','c')
plot(fine.x.psi,fine.y.psi,'k')
plot(fine.x.psi',fine.y.psi','k')
plot(fine.x.rho,fine.y.rho,'ro')
plot(fine.x.u,fine.y.u,'kx')
plot(fine.x.v,fine.y.v,'ms')
 
%to create swan grid then use:
disp(['Created swan grid + bathy files:'])
roms2swan(fine.x.rho,fine.y.rho,fine.h,fine.mask.rho);
