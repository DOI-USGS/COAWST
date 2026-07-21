% plot_duck_xshore
%
% plot results from ducknc xshore
%
% jcw 09Jun2026
%

cd D:\models\COAWST_updates\COAWST_v3.9\Ducknc_xshore_test

% Using 2D spec data at    20231101.00000000
% Peak dir from    90.000000000000000
% Waves hm0 is    2.0020128246191820
% Representative period is    11.674144364983254
% Computing InWave boundary forcing
% Freqs min max are :    2.7777777777777778E-004  0.50000000000000000
% Computing eastern AC forcing
% Mean offshore water depth for AC is    8.9840437021504727

netcdf_load('duck_xshore_his.nc');
figure
plot(x_rho(:,3),-squeeze(bath(:,3,1)),'k','linewidth',2)
hold on
plot(x_rho(:,3),-squeeze(bath(:,3,2:end)))
%
plot(x_rho(:,3),squeeze(zeta(:,3,:)))
zoom on
hold on
plot(x_rho(:,3),squeeze(Hwave(:,3,:)),'--')
title('Ducknc InWave: init and final bathy; 1 mintue water levels (solid) and Hsigs (dashed)')

%cross shore vels
N=size(u,3)
ot=54; %length(ocean_time)
igrid=3 %u
h=squeeze(bath(:,:,ot));
z_u=set_depth(Vtransform, Vstretching, ...
                       theta_s, theta_b, hc, N, ...
                       igrid, h, zeta(:,:,ot));
figure
quiver(repmat(x_u(:,1),1,N),squeeze(z_u(:,3,:)),squeeze(u(:,3,:,ot)),squeeze(u(:,3,:,ot))*0,.01)
zoom on
hold on
plot(x_rho(:,3),-squeeze(bath(:,3,ot)),'k')
plot(x_rho(:,3),squeeze(zeta(:,3,ot)),'b')
ylabel('Depth, m')
xlabel('Cross shore distance, m')
title('Cross shore velocities')


netcdf_load('duck_xshore_avg.nc');
figure
hold on
plot(x_rho(:,3),squeeze(Hwave(:,3,:)))
zoom on
title('Ducknc InWave: 5 min average Hsigs')


