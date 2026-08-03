
%solution to first 2 terms of Oey eq 7 for wet dry slope chan test.
g=9.81;
dx=250;
x=[-125:dx:25125];
r=0.01;  % 0.0025 or 0.01;
depth=10*x/25250;
Hx=10/25250;
D=zeros(2,length(x))+0.10;
%Dsave=zeros(11,length(x));
dt=0.005;
tsave=0;
nnew=1;
nold=2;
mat_time=0;

% 0.5*24*3600=43200 seconds
% dD/dt = -2gDHx/r dD/dx + d/dx(gD^2/r dD/dx)
% dD/dt = -2gDHx/r dD/dx + gD^2/r d/dx(dD/dx) + dD/dx d/dx(gD^2/r)
% dD/dt = -2gDHx/r dD/dx + gD^2/r d/dx(dD/dx) + dD/dx 2Dg/r dD/dx

%save IC
tsave=tsave+1;
Dsave(tsave,:)=D(nnew,:);
mat_time(tsave)=0

if(0)
 for t=dt:dt:43200
  nold=nnew;
  nnew=3-nold;
%  T2=g*D(nold,:).^2/r.*gradient(D(nold,:))/dx;
%  D(nnew,:)=D(nold,:)+dt*(-2*g*D(nold,:)*Hx/r.*gradient(D(nold,:))/dx+gradient(T2)/dx);
%  temp=diff(D(nold,:))/dx;
%  temp(2:end+1)=temp;
%  T2=g*D(nold,1:end).^2/r.*temp;
%  D(nnew,1:end-1)=D(nold,1:end-1)+dt*(-2*g*D(nold,1:end-1)*Hx/r.*diff(D(nold,:))/dx+diff(T2)/dx);

%this is from 2:end-1
   dDdx=(D(nold,3:end)-D(nold,1:end-2))/dx/2;
%   dDdx=dDdx(1:end-1);  %shift to 2:end-1
   d2Ddx2=(D(nold,1:end-2)-2*D(nold,2:end-1)+D(nold,3:end))/dx/dx;
   D(nnew,2:end-1)=D(nold,2:end-1)+dt*(-2*g*D(nold,2:end-1)*Hx/r.*dDdx+g.*D(nold,2:end-1).^2/r.*d2Ddx2+ ...
                   dDdx.*dDdx*2.*D(nold,2:end-1)*g/r);
   D(nnew,end)=10*sin(2*pi/(12*3600)*t/2)+0.10;
   if (((1-mod(t,4320))==1)&(t>0))
     tsave=tsave+1;
     Dsave(tsave,:)=D(nnew,:);
     mat_time(tsave)=t
   end
 end
end

%save wetdry_slope_chan0025.mat

%plotting
if 1
  load e:\data\matlab\test_cases\wetdry\wetdry_slope_chan0025.mat
%  ncload c:/work/models/COAWST/wetdry_slope_chan0025_his.nc
%  ncload d:\data\test_cases\wetdry\ocean_wetdry_slope_chan0025_his.nc
 % netcdf_load e:\data\models\COAWST\wetdry_slope_chan0025_his.nc
  netcdf_load e:\data\matlab\test_cases\wetdry\wetdry_slope_chan0025_his.nc
  figure
  subplot(211)
  hold on
  clear rmse1
  for mm=[2 4 6 7 9 11]
%   plot(x_rho(:,2)/1000,squeeze(zeta(:,3,mm))+h(:,4),'k--')
%   plot(x/1000,Dsave(mm,:)','k')
    plot(x_rho(:,2)/1000,squeeze(zeta(:,3,mm)),'k--')
    plot(x/1000,Dsave(mm,:)'-h(:,4),'k')
    a=squeeze(zeta(:,3,mm));
    b=Dsave(mm,:)'-h(:,4);
    rmse1(mm)=sqrt(sum((a-b).^2)/length(a));
  end
% set(gca,'xlim',[0 25],'ylim',[-1 11])
  plot(x/1000,-h(:,4),'r')
  set(gca,'xlim',[0 25],'ylim',[-11 1])
  title('r=0.0025')

%  clear
  load e:\data\matlab\test_cases\wetdry\wetdry_slope_chan01.mat
%  ncload c:/work/models/COAWST/wetdry_slope_chan01_his.nc
%  ncload d:\data\test_cases\wetdry\ocean_wetdry_slope_chan01_his.nc
%  netcdf_load e:\data\models\COAWST\wetdry_slope_chan01_his.nc
  netcdf_load e:\data\matlab\test_cases\wetdry\wetdry_slope_chan01_his.nc
  subplot(212)
  hold on
  clear rmse2
  for mm=[2 4 6 7 9 11]
%   plot(x_rho(:,2)/1000,squeeze(zeta(:,3,mm))+h(:,4),'k--')
%   plot(x/1000,Dsave(mm,:)','k')
    plot(x_rho(:,2)/1000,squeeze(zeta(:,3,mm)),'k--')
    plot(x/1000,Dsave(mm,:)'-h(:,4),'k')
    a=squeeze(zeta(:,3,mm));
    b=Dsave(mm,:)'-h(:,4);
    rmse2(mm)=sqrt(sum((a-b).^2)/length(a));
  end
% set(gca,'xlim',[0 25],'ylim',[-1 11])
  plot(x/1000,-h(:,4),'r')
  set(gca,'xlim',[0 25],'ylim',[-11 1])
  title('r=0.01')
end


