clear all
close all

hisfile='ocean_his.nc';
h=squeeze(ncread(hisfile,'h')); 
X=squeeze(ncread(hisfile,'x_rho')); 
Y=squeeze(ncread(hisfile,'y_rho')); 
time=ncread(hisfile,'ocean_time');
nt=length(time);

count=0;
for tt=1:1:nt
  count=count+1;
% mask=squeeze(ncread(hisfile,'wetdry_mask_rho',[1 1 tt],[Inf Inf 1]));
  mask=squeeze(ncread(hisfile,'mask_rho',[1 1 tt],[Inf Inf 1]));
  zeta=squeeze(ncread(hisfile,'zeta',[1 1 tt],[Inf Inf 1]));
  Hwave=squeeze(ncread(hisfile,'Hwave',[1 1 tt],[Inf Inf 1]));
  dum=find(mask==0);
  zeta(dum)=nan;
  Hwave(dum)=nan;

  figure(2)
  set(gcf,'renderer','painters')
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
  fnames=strcat('Hwave_zeta_',num2str(count),'_0_1_DEAN');
  print('-dpng','-r200',fnames)
end
% 

figure
zetaz=squeeze(ncread(hisfile,'zeta')); 
zetaz=squeeze(zetaz(:,150,:)); 
wetdry=squeeze(ncread(hisfile,'wetdry_mask_rho',[1 150 1],[Inf 1 Inf]));
%zetaz(wetdry==0)=nan;
amp=max(zetaz(:,:).')-mean(zetaz(:,:).');
amp1=min(zetaz(:,:).');
plot(X(:,150),amp,'Linewidth',1)
hold on
plot(X(:,150),mean(zetaz.'),'r','Linewidth',1)
xlabel('Cross shore direction')
ylabel('\eta amplitude (m)')
legend('Max water level -setup','Setup or setdown')
axis([0 970 0 1])
set(gcf,'PaperPosition',([0.1 0.1 6 3]))
fnames=strcat('amp_zeta_0_1_DEAN');
print('-dpng','-r200',fnames)

