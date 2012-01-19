clear all
close all

hisfile='P:\INWAVE_020811\Projects\Inwave_tests\Inwave_shoreface\ocean_his.nc';
%hisfile='P:\INWAVE_020811\Projects\Inwave_tests\Inwave_shoreface\ocean_his_0_0015.nc';
nc=netcdf(hisfile);
h=squeeze(nc{'h'}(:,:)); 
X=squeeze(nc{'x_rho'}(:,:)); 
Y=squeeze(nc{'y_rho'}(:,:)); 
zetaz=squeeze(nc{'zeta'}(:,300,20)); 

time=squeeze(nc{'ocean_time'}(:,:));
nt=length(time);

count=0;
for tt=1:1:nt
count=count+1;
mask=squeeze(nc{'wetdry_mask_rho'}(tt,:,:));
zeta=squeeze(nc{'zeta'}(tt,:,:));
Hwave=squeeze(nc{'Hwave'}(tt,:,:));
dum=find(mask==0);
zeta(dum)=nan;
Hwave(dum)=nan;

figure(2)
subplot(1,2,1)
pcolor(X,Y,Hwave)
shading flat
colorbar 'vertical'
caxis([0 3.5])
axis([0 1070 0 6000])
title('Hs')
subplot(1,2,2)
pcolor(X,Y,zeta)
shading flat
colorbar 'vertical'
title('\eta')
caxis([-0.4 1.2])
axis([0 1070 0 6000])

set(gcf,'PaperPosition',([0.1 0.1 6 7]))
fnames=strcat('P:\INWAVE_020811\Projects\Inwave_tests\Inwave_shoreface\Hwave_zeta_',num2str(count),'_0_1_DEAN');
print('-dpng','-r200',fnames)
end
% 

figure
hisfile='P:\INWAVE_020811\Projects\Inwave_tests\Inwave_shoreface\ocean_his.nc';
nc=netcdf(hisfile);
zetaz=squeeze(nc{'zeta'}(:,150,:)); 
wetdry=squeeze(nc{'wetdry_mask_rho'}(:,150,:));
zetaz(wetdry==0)=nan;
amp=max(zetaz(100:300,:))-mean(zetaz(100:300,:));
amp1=min(zetaz(100:300,:));
plot(X(150,:),amp,'Linewidth',1)
hold on
plot(X(150,:),mean(zetaz),'r','Linewidth',1)
xlabel('Cross shore direction')
ylabel('\eta amplitude (m)')
legend('Max water level -setup','Setup or setdown')
axis([0 970 0 1])
set(gcf,'PaperPosition',([0.1 0.1 6 3]))
fnames=strcat('P:\INWAVE_020811\Projects\Inwave_tests\Inwave_shoreface\amp_zeta_0_1_DEAN');
print('-dpng','-r200',fnames)

% kk(:,:)=amp(1,:,:);
% pcolor(kk)
% shading flat

% shading flat
% colorbar 'vertical'
