clear all
close all

grdfile='N:\INWAVE_011110\Projects\Inwave_test_mai\beach\InWave_grd.nc';
nc=netcdf(grdfile);
h=squeeze(nc{'h'}(10,:)); 
prof=squeeze(nc{'h'}(:,:)); 
mask=squeeze(nc{'mask_rho'}(:,:)); 
dum=find(mask==0);
prof(dum)=nan;
close (nc)

[X,Y]=meshgrid(0:100:60*100, 0:100:19*100);

figure(2)
pcolor(X,Y,prof)
shading flat
xlabel('Distance, m')
ylabel('Distance, m')
colorbar 'vertical'

set(gcf,'PaperPosition',([0.1 0.1 5 2]))
fnames='N:\INWAVE_011110\Projects\Inwave_test_mai\beach\bathy';
print('-dpng','-r200',fnames)

bndfile='N:\INWAVE_011110\Projects\Inwave_test_mai\beach\InWave_bnd.nc';
nc=netcdf(bndfile);
time_bound=squeeze(nc{'energy_time'}(:)); 
AC_bound=squeeze(nc{'AC_west'}(:,1,10)); 
close (nc)

figure(3)
plot(time_bound,AC_bound,'LineWidth',2)
xlabel('Time, sec')
ylabel('Action Balance')
set(gcf,'PaperPosition',([0.1 0.1 5 2]))
fnames='N:\INWAVE_011110\Projects\Inwave_test_mai\beach\boundary';
print('-dpng','-r200',fnames)

hisfile='N:\INWAVE_011110\Projects\Inwave_test_mai\beach\ocean_his.nc';
nc=netcdf(hisfile);
AC_nb=squeeze(nc{'AC'}(:,16,:,:)); 
E_nb=(2*pi/10).*max(AC_nb);
H_nb=(8.*abs(E_nb)/(1028.2845*9.81)).^0.5;
close (nc)


hisfile='N:\INWAVE_011110\Projects\Inwave_test_mai\beach_bl\ocean_his.nc';
nc=netcdf(hisfile);
AC_bl=squeeze(nc{'AC'}(:,16,:,:)); 
E_bl=(2*pi/10).*max(AC_bl);
H_bl=(8.*abs(E_bl)/(1028.2845*9.81)).^0.5;
close (nc)

hisfile='N:\INWAVE_011110\Projects\Inwave_test_mai\beach_br\ocean_his.nc';
nc=netcdf(hisfile);
AC_br=squeeze(nc{'AC'}(:,16,:,:)); 
E_br=(2*pi/10).*max(AC_br);
H_br=(8.*abs(E_br)/(1028.2845*9.81)).^0.5;
close (nc)

figure (1)

subplot(3,2,1)
AC_nb1(:,:)=AC_nb(:,10,:);
AC_nb1(AC_nb1>1000000)=nan;
pcolor(AC_nb1)
shading flat
xlabel('Cell number in the xi direction')
ylabel('Time, minutes')
colorbar 'vertical'
title('Action balance, no dissipation')
subplot(3,2,2)
H_nb1(:,:)=H_nb(1,:,:);
pcolor(H_nb1.*mask)
shading flat
xlabel('Cell number in the xi direction')
ylabel('Cell number in the etai direction')
colorbar 'vertical'
title('Maximum Wave Height')
% 
subplot(3,2,3)
AC_br1(:,:)=AC_br(:,10,:);
AC_br1(AC_br1>1000000)=nan;
pcolor(AC_br1)
shading flat
xlabel('Cell number in the xi direction')
ylabel('Time, minutes')
colorbar 'vertical'
title('Action balance, Roelvink')
subplot(3,2,4)
H_br1(:,:)=H_br(1,:,:);
pcolor(H_br1.*mask)
shading flat
xlabel('Cell number in the xi direction')
ylabel('Cell number in the etai direction')
colorbar 'vertical'
title('Maximum Wave Height')

subplot(3,2,5)
AC_bl1(:,:)=AC_bl(:,10,:);
AC_bl1(AC_bl1>1000000)=nan;
pcolor(AC_bl1)
shading flat
xlabel('Cell number in the xi direction')
ylabel('Time, minutes')
colorbar 'vertical'
title('Action balance, H=\gamma * h')
subplot(3,2,6)
H_bl1(:,:)=H_bl(1,:,:);
pcolor(H_bl1.*mask)
shading flat
xlabel('Cell number in the xi direction')
ylabel('Cell number in the etai direction')
colorbar 'vertical'
title('Maximum Wave Height')

set(gcf,'PaperPosition',([0.1 0.1 8 8]))
fnames='N:\INWAVE_011110\Projects\Inwave_test_mai\beach\solution';
print('-dpng','-r200',fnames)

figure(5)

gamma_nb=H_nb1./prof;
gamma_bl=H_bl1./prof;
gamma_br=H_br1./prof;
    

subplot(2,1,1)
pcolor(gamma_bl)
shading flat
xlabel('Cell number in the xi direction')
ylabel('Time, minutes')
colorbar 'vertical'
title('H/h, H=\gamma * h')
subplot(2,1,2)
pcolor(gamma_br)
shading flat
xlabel('Cell number in the xi direction')
ylabel('Time, minutes')
colorbar 'vertical'
title('H/h, Roelvink')

set(gcf,'PaperPosition',([0.1 0.1 8 8]))
fnames='F:\MAITANE\PROJECTS\INWAVE\COAWST_MAI\Projects\Inwave_tests\beach\gamma';
print('-dpng','-r200',fnames)

figure(6)
% plot(X(10,:),h,'--')
hold on
% plot(X(10,:),gamma_nb(10,:))
plot(X(10,:),gamma_br(10,:),'g','LineWidth',2)
plot(X(10,:),gamma_bl(10,:),'r','LineWidth',2)
ylabel('\gamma')
xlabel('Distance in the xi direction')
max_H(1,:)=H_nb(1,20,:);
% plot(X(10,:),max_H,'c')
legend('Run 1','Run 2','Run 3')


set(gcf,'PaperPosition',([0.1 0.1 4 4]))
fnames='F:\MAITANE\PROJECTS\INWAVE\COAWST_MAI\Projects\Inwave_tests\beach\profilr';
print('-dpng','-r200',fnames)

