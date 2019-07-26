% mfile to plot results from Rip_current tets case
% uses native matlab netcdf interface
% jcwarner 05Spet2013
%

%enter file name
ncfile='ocean_rip_current_his.nc';

%get some vars
ubar=ncread(ncfile,'ubar');
vbar=ncread(ncfile,'vbar');
x_rho=ncread(ncfile,'x_rho');
y_rho=ncread(ncfile,'y_rho');
h=ncread(ncfile,'h');

%avg data to similar points
us=squeeze(ubar(:,:,end));
vs=squeeze(vbar(:,:,end));
us=0.5*(us(1:end-1,2:end-1)+us(2:end,2:end-1));
vs=0.5*(vs(2:end-1,1:end-1)+vs(2:end-1,2:end));
xs=x_rho(2:1:end-1,2:1:end-1);
ys=y_rho(2:1:end-1,2:1:end-1);

%plot it out
figure
pcolorjw(x_rho,y_rho,h);
colorbar
hold on
inc=3;
quiver(xs(1:inc:end,1:inc:end),ys(1:inc:end,1:inc:end), ...
       us(1:inc:end,1:inc:end),vs(1:inc:end,1:inc:end),3,'w')
xlabel('Distance (m)')
ylabel('Distance (m)')
title('Bathymetry (color) and depth-avg velocity (arrows)')

%print -dpng rip_current_vels.png


