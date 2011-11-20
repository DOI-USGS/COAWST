clear all
close all

path_1='F:\MAITANE\PROJECTS\INWAVE\COAWST_MAI\Projects\Inwave_tests\bicro\';

archivo=strcat(path_1,'ocean_his.nc');
nc=netcdf(archivo); 
h=squeeze(nc{'h'}(:,:)); 
dum=find(h<0);

for t=1:95
zeta=squeeze(nc{'zeta'}(t,:,:));  
zeta(dum)=nan;
figure
pcolor(zeta)
shading flat
caxis([-0.005 0.005])
pause
clear zeta
close (1)
end

for i=1:1:131
zeta1=squeeze(nc{'zeta'}(:,10,i)); 
plot(zeta1)
axis([0 100 -0.006 0.006])
pause
end