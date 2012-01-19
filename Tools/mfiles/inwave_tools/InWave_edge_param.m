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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%        GENERAL PARAMETERS: THESE NEED TO BE DEFINED ALWAYS       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lm= 140;               % number of total rho points in the xi direction
Mm= 500;               % number of total rho points in the eta direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                GRID AND BATHYMETRY DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath='P:\INWAVE_020811\Projects\Inwave_tests\Edge\';

if (make_InWave_grd)
    
  grd_file=strcat(filepath,'InWave_grd.nc');  % name of the grid file
    
  % Grid characteristics

  dx=5;                % grid cell size in xi direction
  dy=5;                % grid cell size in eta direction 

  % enter x and y coordinates of rho points

  x=[-dx/2:dx:dx*(Lm-1)];
  y=[-dy/2:dy:dy*(Mm-1)];

  [x,y]=meshgrid(x,y);
  
%   repmat(x,length(y),1);
%   y=repmat(y',length(y),1);

  % Bathymetry characteristics
  
  depth0= 18;            % water depth in the study domain (m)
  angle= 0.03;

  % set depth 
  
  depth=zeros(size(x))+depth0;
  
  k=2*pi/500;
  Amp=0.05;

  for jj=1:Mm
      for ii=1:Lm
          depth(jj,ii)=18.0-(140-ii)*dx*angle+Amp*cos(k*jj*dy);
      end
  end

  for jj=1:100
      for ii=1:Lm
          depth(jj,ii)=18.0-(140-ii)*dx*angle+Amp;
      end
  end

  for jj=401:500
      for ii=1:Lm
          depth(jj,ii)=18.0-(140-ii)*dx*angle+Amp;
      end
  end

  depth(depth<0.15)=0.15;
  pcolor(depth)
  shading flat
  colorbar 'vertical'

  % set grid angle
  roms_angle=zeros(size(depth));

  % set masking
  mask_rho=ones(size(depth));
  dum=find(depth<=0.15);
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

  Ac=ones(Nbins,Mm,Lm).*0;
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
  
  dt= 1.000;             % time increment for the boundaries (seconds)
  drtn= 3600;            % total duration of the simulation (seconds)

  time=[0:dt:drtn];
  
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
    
    Hrms=0.75; 
fIG=0.009;
Tp=10;

a2=sqrt((Hrms/2).^2/5);
a1=2*a2;
T1=Tp;
f1=1/T1;
f2=f1-fIG;
T2=1./f2;
Tm=(T1+T2)./2;
TA=Tm;

eta=a1*cos(2*pi*f1*time)+a2*cos(2*pi*f2*time);

en=hilbert(eta);

E=1/8*1025*9.81*(2.*abs(hilbert(eta))).^2;
Ac1=E./(2*pi/Tm);

plot(time,eta)
hold on
plot(time,abs(hilbert(eta)),'g')

    for i=1:Mm
    Ac_east(:,1,i)=Ac1(:);
    end
    
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



