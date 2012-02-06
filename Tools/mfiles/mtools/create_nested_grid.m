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
% Version jcwarner Aug 3, 2007
% Updated clea denamiel Feb 5, 2012

%%%%%%%%%%%%% START OF USER SECTION %%%%%%%%%%%%%%%%

%1) ENTER NAME OF EXISITNG COARSE GRID FILE
%ncfile_coarse='refined_chan_grid.nc';
ncfile_coarse='inlet_test_grid.nc';

%2) ENTER NAME OF NEW FINE GRID FILE
%ncfile_fine='refined_chan_grid_ref5.nc';
ncfile_fine='inlet_test_grid_ref5_test.nc';

%3) ENTER THE START AND END INDICES OF THE PSI POINTS OF THE COARSE GRID
%   THAT IDENTIFY THE OUTER BOUNDS OF THE FINE GRID
%Istr=30; Iend=50; Jstr=2; Jend=5;  % test_chan_refined
Istr=24; Iend=54; Jstr=40; Jend=56;  % inlet_test_refined

%4) ENTER SCALE FACTOR FOR INCREASED RESOLUTION (use 3 or 5)
scale=5;

%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

%set some random variables here
msize1=5;      %marker size for plots

%get the coarser grid
netcdf_load(ncfile_coarse)

%number of rho points in coarse grid section to have increased resolution
Ilen_c=Iend-Istr;
Jlen_c=Jend-Jstr;

%Now set dimensions for fine grid.
%xi direction (left-right)
LP = Ilen_c*scale+7;   %rho dimension
L = LP-1;              %psi dimension

%eta direction (up-down)
MP = Jlen_c*scale+7;   % rho dimension
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

%set some arrays for interpolation here, psi points
% _c for coarse, _f for fine
offset=ceil(4/scale);  % need 4 points to the left, 3 to the right
numx_c=((Iend+1)-(Istr-offset))+1;
numy_c=((Jend+1)-(Jstr-offset))+1;
numx_f=scale*(numx_c-1)+1;
numy_f=scale*(numy_c-1)+1;

xc_c=[1:numx_c];
xc_c=repmat(xc_c,numy_c,1);
yc_c=[1:numy_c];
yc_c=repmat(yc_c',1,numx_c);

xc_f=[1:1/scale:numx_c];
xc_f=repmat(xc_f,numy_f,1);
yc_f=[1:1/scale:numy_c];
yc_f=repmat(yc_f',1,numx_f);

%%%%%%%%%  Lon Lat SPACE  %%%%%%%%%%%
%establish a full grid of psi points, including an extra row and column
x_full_grid=interp2(xc_c, yc_c, lon_psi(Jstr-offset:Jend+1,Istr-offset:Iend+1), xc_f, yc_f);
y_full_grid=interp2(xc_c, yc_c, lat_psi(Jstr-offset:Jend+1,Istr-offset:Iend+1), xc_f, yc_f);

%set the psi points
fine.lon.psi=x_full_grid(offset+2:end-(scale-3+1),offset+2:end-(scale-3+1));
fine.lat.psi=y_full_grid(offset+2:end-(scale-3+1),offset+2:end-(scale-3+1));

%set the rho points
ztemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lon.rho=ztemp(2:2:end-1,2:2:end-1);
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lat.rho=ztemp(2:2:end-1,2:2:end-1);

%set the u points
ztemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lon.u=ztemp(2:2:end-1,3:2:end-2);
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lat.u=ztemp(2:2:end-1,3:2:end-2);

%set the v points
ztemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lon.v=ztemp(3:2:end-2,2:2:end-1);
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lat.v=ztemp(3:2:end-2,2:2:end-1);

% compute dx and dy at u and v points for future computation of
% pm, pn, dmde, dndx
xtemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
ytemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
dx = earthdist(xtemp(:, 2:end), ytemp(:, 2:end), xtemp(:, 1:end-1), ytemp(:, 1:end-1));
dy = earthdist(xtemp(2:end, :), ytemp(2:end, :), xtemp(1:end-1, :), ytemp(1:end-1, :));
sx = 0.5*(dx(1:end-1, :) + dx(2:end, :));
sy = 0.5*(dy(:, 1:end-1) + dy(:, 2:end));
sx= sx(1:2:end-1,1:2:end-1);
sy= sy(2:2:end,2:2:end);
fine.pm = 1 ./ sx;
fine.pn = 1 ./ sy;
fine.pm(isinf( fine.pm))=0.999e-3;
fine.pn(isinf( fine.pn))=0.999e-3;
fine.pm(isnan( fine.pm))=0.999e-3;
fine.pn(isnan( fine.pn))=0.999e-3;
fine.dmde = zeros(size(fine.pm));
fine.dndx = zeros(size(fine.pn));
fine.dmde(2:end-1, :) = 0.5*(1./fine.pm(3:end, :) - 1./fine.pm(1:end-2, :));
fine.dndx(:, 2:end-1) = 0.5*(1./fine.pn(:, 3:end) - 1./fine.pn(:, 1:end-2));
fine.dmde(isinf(fine.dmde))=0;
fine.dndx(isinf(fine.dndx))=0;
fine.dmde(isnan(fine.dmde))=0;
fine.dndx(isnan(fine.dndx))=0;
% Grid-cell orientation, degrees counter-clockwise from
%  east, presently based on flat-earth approximation.
RCF=180/pi;
dlon_temp = diff(xtemp.').';
dlat_temp = diff(ytemp.').';
clat = cos(ytemp / RCF);
clat(:, end) = [];
ang = atan2(dlat_temp, dlon_temp .* clat);
temp = 0.5*(ang(1:end-1, :) + ang(2:end, :));
temp(isnan(temp))=0; % doesn't like NaN
fine.angle= temp(1:2:end,1:2:end);

%%%%%%%%%  X Y SPACE  %%%%%%%%%%%
%establish a full grid of psi points, including an extra row and column
x_full_grid=interp2(xc_c, yc_c, x_psi(Jstr-offset:Jend+1,Istr-offset:Iend+1), xc_f, yc_f);
y_full_grid=interp2(xc_c, yc_c, y_psi(Jstr-offset:Jend+1,Istr-offset:Iend+1), xc_f, yc_f);

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
fine.x.u=ztemp(2:2:end-1,3:2:end-2);
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.y.u=ztemp(2:2:end-1,3:2:end-2);

%set the v points
ztemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.x.v=ztemp(3:2:end-2,2:2:end-1);
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.y.v=ztemp(3:2:end-2,2:2:end-1);


%%%%%% GRID Metrics  %%%%%%%%%

[XI,YI]=meshgrid([min(x_psi(:)):min(1./pm(:))/2:max(x_psi(:))], ...
                 [min(y_psi(:)):min(1./pn(:))/2:max(y_psi(:))]);
ZI=griddata(x_rho,y_rho,h,XI,YI);
%ZI=interp2(x_rho,y_rho,h,XI,YI);
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
%%%%%%%%%%%%%%%%%%% finished creating grid %%%%%%%%%%%%%%%%

%call routine to create netcdf file
save('test.mat','-struct','fine')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot new grid inside old grid
figure(1)
plot(x_psi,y_psi,'c')
hold on
plot(x_psi',y_psi','c')
plot(fine.x.psi,fine.y.psi,'k')
plot(fine.x.psi',fine.y.psi','k')
plot(fine.x.rho,fine.y.rho,'ro')
plot(fine.x.u,fine.y.u,'kx')
plot(fine.x.v,fine.y.v,'ms')
plot(x_rho,y_rho,'bo')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bathymetry
figure(2)
subplot(2,1,1);
pcolor(x_rho,y_rho,h)
shading flat
colorbar
title('Parent grid - bathymetry (m)');
subplot(2,1,2);
pcolor(fine.x.rho,fine.y.rho,fine.h)
shading flat
colorbar
title('Child grid - bathymetry (m)');
% PM
figure(3)
subplot(2,1,1);
pcolor(x_rho,y_rho,pm)
shading flat
colorbar
title('Parent grid - pm');
subplot(2,1,2);
pcolor(fine.x.rho,fine.y.rho,fine.pm)
shading flat
colorbar
title('Child grid - pm');
% PN
figure(4)
subplot(2,1,1);
pcolor(x_rho,y_rho,pn)
shading flat
colorbar
title('Parent grid - pn');
subplot(2,1,2);
pcolor(fine.x.rho,fine.y.rho,fine.pn)
shading flat
colorbar
title('Child grid - pn');
% dndx
figure(5)
subplot(2,1,1);
pcolor(x_rho,y_rho,dndx)
shading flat
colorbar
title('Parent grid - dndx');
subplot(2,1,2);
pcolor(fine.x.rho,fine.y.rho,fine.dndx)
shading flat
colorbar
title('Child grid - dndx');
% dmde
figure(6)
subplot(2,1,1);
pcolor(x_rho,y_rho,dmde)
shading flat
colorbar
title('Parent grid - dmde');
subplot(2,1,2);
pcolor(fine.x.rho,fine.y.rho,fine.dmde)
shading flat
colorbar
title('Child grid - dmde');
% angle
figure(7)
subplot(2,1,1);
pcolor(x_rho,y_rho,angle)
shading flat
colorbar
title('Parent grid - angle (rad)');
subplot(2,1,2);
pcolor(fine.x.rho,fine.y.rho,fine.angle)
shading flat
colorbar
title('Child grid - angle (rad)');

%to create swan grid then use:
disp(['Created swan grid + bathy files:'])
roms2swan(fine.x.rho,fine.y.rho,fine.h,fine.mask.rho);
