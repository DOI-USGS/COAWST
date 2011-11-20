% THIS FILE CONTAINS THE DEFINITIONS OF THE PARAMETERS NEEDED TO CREATE THE 
% INPUT FILES FOR THE INWAVE MODEL:
%
% InWave_grd.nc
% InWave_ini.nc
% InWave_bnd.nc
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                DEFINE WHICH FILES TO GENERATE                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_InWave_grd=1;
make_InWave_ini=1;
make_InWave_bnd=1;


%%% LOAD DUCK85 DATA

%% LOAD THE SEA SURFACE ELEVATION DATA AND GET THE WAVE ENVELOPE

eta=load('F:\MAITANE\PROJECTS\INWAVE\DUCK_85\d0915h5e.asc');
t=[0:0.5:(4800-1).*0.5];
time=t;

figure(2)
plot(t,eta)

y=eta;
df=1/(t(2)-t(1));
a=fft(y);
a=a';
f1=[linspace(0,df/2,length(a)/2+1)];
amp=2.*abs(a(1:length(f1)))/length(t);
phase=angle(a(1:length(f1)))';

clear y

% Filter frequencies higher than 0.1 Hz
amp1=amp;
dum=find(f1>0.25);
amp1(dum)=0;

nf=length(f1);
nt=length(t);
eta1=zeros(1,nt);

for ii=1:nt
  for jj=1:nf
    eta1(1,ii)=eta1(1,ii)+amp1(1,jj)*cos(-2*pi*f1(1,jj)*(ii-1)*0.5+phase(jj,1));
  end
end

figure(5)
plot(t,eta)
hold on
plot(time,eta1,'c')

x_hil = hilbert(eta1);
amp_hil=(real(x_hil).^2+imag(x_hil).^2).^0.5;

plot(t,amp_hil,'g')
xlabel('Time, sec')
ylabel('\eta, m')


sea=[t',eta1'];
count=0;
for i=1:length(sea(:,2))-1
    if(sea(i+1,2)<=0 && sea(i,2)>0)
    count=count+1;
    pos(count)=i;
    end
end
nwaves=count;
count=0;
for i=1:nwaves-1
    minh=min(sea(pos(i):pos(i+1),2));
    maxh=max(sea(pos(i):pos(i+1),2));
    H(i)=abs(maxh)+abs(minh);
    T(i)=t(pos(i+1))-t(pos(i));
    time_w(i)=(t(pos(i+1))+t(pos(i)))/2;
end
figure(11)
hold on
plot(t,eta1,'k')
plot(time_w,H,'.r') 
plot(t,amp_hil,'g')
HH=interp1(time_w,H,t,'pchip');
hold on
plot(t,HH./2,'b')
 
plot(t,-amp_hil,'g')
plot(t,-HH./2,'b')
axis([0 300 -0.6 0.6])

legend('Smoothed \eta','Individual upcrossing H','Hilbert','Interpolated H/2 envelope',4)
axis([0 300 -0.7 0.7])
str1=strcat('F:\MAITANE\PROJECTS\INWAVE\DUCK_85\input');
pp=[0.1 0.1 10 4];
set(gcf,'paperposition',pp)
set(gcf,'position',pp)
print('-dpng','-r600',str1)


%% LOAD THE PROFILE DATA

load('F:\MAITANE\PROJECTS\INWAVE\DUCK_85\p091030.dat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%        GENERAL PARAMETERS: THESE NEED TO BE DEFINED ALWAYS       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lm= 71;                % number of total rho points in the xi direction
Mm= 20;        %100    % number of total rho points in the eta direction
TA= 11;                % representative absolute wave period (sec)

%
amp_HH=(1/8*1025*9.81*(HH).^2)./(2*pi/TA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                GRID AND BATHYMETRY DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath='..\..\Projects\Inwave_test_mai\Duck_85\';


if (make_InWave_grd)
    
  grd_file=strcat(filepath,'InWave_grd.nc');  % name of the grid file
    
  % Grid characteristics

  dx=5;                % grid cell size in xi direction
  dy=5;                % grid cell size in eta direction 

  % enter x and y coordinates of rho points

  x1=[0:dx:dx*(Lm-1)];
  y1=[0:dy:dy*(Mm-1)];

  x=repmat(x1,length(y1),1);
  y=repmat(y1',1,length(x1));

  % Bathymetry characteristics
  depth0= p091030(:,2)';            % water depth in the study domain (m)

  % set depth 
  depth=repmat(depth0,length(y1),1);
  dum=find(depth<=0.2);
  depth(dum)=0.2;
  
  % set grid angle
  roms_angle=zeros(size(depth));

  % set masking
  mask_rho=ones(size(depth));
  mask_rho(dum)=0;

  % set coriolis f
  f=zeros(size(depth))+4.988e-5; %20N

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 INITIAL CONDITION DEFINITION                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_ini || make_InWave_bnd )  

  Nbins= 20;                               % number of directional bins considered in the simulation
  Bindirs = [0:360/Nbins:360-360/Nbins];   % center angles of the directional bins.
  pd=ones(size(Bindirs)).*360./(Nbins);
  Bindirs_c = [0-(360/Nbins)/2:360/Nbins:360-360/Nbins+(360/Nbins)/2];

end  

if (make_InWave_ini)  
    
  ini_file=strcat(filepath,'InWave_ini.nc');  % name of the initial file

  Ac=ones(Nbins,Mm,Lm).*0.0001;
  Cx=ones(Nbins,Mm,Lm-1).*0;
  Cy=ones(Nbins,Mm-1,Lm).*0;
  Ct=ones(Nbins+1,Mm,Lm).*0;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 BOUNDARY CONDITION DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_bnd)

  bnd_file=strcat(filepath,'InWave_bnd.nc');  % name of the boundary file

  % Duration of the simulation and time increment for the boundaries
   
  dt= t(1,2)-t(1,1);       % time increment for the boundaries (seconds)
  drtn= t(1,end);          % total duration of the simulation (seconds)

 
  
  % Specify by 1 the open boundary: S E N W
  obc=[0 1 0 0];

%    Nbins_bnd= 20;         % number of directional bins with energy at the boundaries
%    dir_bnd= Bindirs;      % center angle (degrees of the bin containing energy at the boudaries
 
 % Specify number of directions at the boundaries (we have to specify at
 % least 2)
 
  Nbins_bnd= 1;         % number of directional bins with energy at the boundaries
  dir_bnd= 90;         % center angle (degrees of the bin containing energy at the boudaries)

  if sum(ismember(dir_bnd,Bindirs)) ~= Nbins_bnd; 
  bin_error=1;
  else
  bin_error=0;
  end      

  if obc(1)==1
    Ac_north=zeros(length(time),Nbins_bnd,Lm);
  end
  if obc(2)==1
    Ac_east=zeros(length(time),Nbins_bnd,Mm);
    amp_HH1=repmat(amp_HH,Mm,1)';
    Ac_east(:,1,:)=amp_HH1(:,:);
  end
  if obc(3)==1
    Ac_south=zeros(length(time),Nbins_bnd,Lm);
  end
  if obc(4)==1
    Ac_west=zeros(length(time),Nbins_bnd,Mm);   
  end

  if obc(1)==1
    TA_west=TA;
  end
  if obc(2)==1
    TA_east=TA;
  end
  if obc(3)==1
    TA_south=TA;
  end
  if obc(4)==1
    TA_north=TA;
  end
end

if (make_InWave_grd)
bry_file=strcat(filepath,'InWave_bry.nc'); 

y1=eta(:);
df=1/(t(2)-t(1));
a=fft(y1);
a=a';
f1=[linspace(0,df/2,length(a)/2+1)];
amp=2.*abs(a(1:length(f1)))/length(t);
phase=angle(a(1:length(f1)))';

% Filter frequencies higher than 0.1 Hz
amp1=amp;
dum=find(f1>0.05);
amp1(dum)=0;

nf=length(f1);
nt=length(t);
zeta=zeros(1,nt);

for ii=1:nt
  for jj=1:nf
    zeta(1,ii)=zeta(1,ii)+amp1(1,jj)*cos(-2*pi*f1(1,jj)*(ii-1)*0.5+phase(jj,1));
  end
end

ubar=load('F:\MAITANE\PROJECTS\INWAVE\DUCK_85\d0915h5u.asc');
y1=ubar(:);
df=1/(t(2)-t(1));
a=fft(y1);
a=a';
f1=[linspace(0,df/2,length(a)/2+1)];
amp=2.*abs(a(1:length(f1)))/length(t);
phase=angle(a(1:length(f1)))';

% Filter frequencies higher than 0.1 Hz
amp1=amp;
dum=find(f1>0.05);
amp1(dum)=0;

nf=length(f1);
nt=length(t);
ubar=zeros(1,nt);

for ii=1:nt
  for jj=1:nf
    ubar(1,ii)=ubar(1,ii)+amp1(1,jj)*cos(-2*pi*f1(1,jj)*(ii-1)*0.5+phase(jj,1));
  end
end

end


