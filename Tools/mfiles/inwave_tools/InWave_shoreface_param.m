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

Lm= 122;                % number of total rho points in the xi direction
Mm= 302;                % number of total rho points in the eta direction
TA= 8.3;                % representative absolute wave period (sec)
theta=261;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                GRID AND BATHYMETRY DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath='P:\INWAVE_020811\Projects\Inwave_tests\Inwave_shoreface\';

if (make_InWave_grd)
    
  grd_file=strcat(filepath,'InWave_grd.nc');  % name of the grid file
    
  % Grid characteristics
  x1=[1:Lm];
  y1=[1:Mm];
  dx=ones(size(x1));
  dy=ones(size(y1));
  dx(1:30)=20;
  dx(31:61)=20-((31:61)-30)./30*18;
  dx(62:end)=2;
  dy=20.*dy;
%   dx=20;                % grid cell size in xi direction
%   dy=20;                % grid cell size in eta direction 
   count=0;
   for ii=1:length(dx)
       if ii>1
           X(ii)= count+dx(ii);
           count= X(ii);
       else
           X(ii)= count;
           count= X(ii);
       end
   end
  % enter x and y coordinates of rho points

%   X=[0:dx:dx*(Lm-1)];
  Y=[0:dy:dy*(Mm-1)];

  [x,y]=meshgrid(X,Y);
  [dx,dy]=meshgrid(dx,dy);
  %x=repmat(x,length(y),1);
  %y=repmat(y',1,length(x));

  % Bathymetry characteristics
  depth0= 11.75;            % water depth in the study domain (m)
  slope=0.0125;
%  slope=0.125;

  % set depth 
%   depth=zeros(size(x))+depth0;
%   depth=depth0-slope.*x;
  depth=0.125*(-(x-max(max(x)))).^(2/3)-3;
  depth(depth<-2)=-2;
  %depth0-slope.*x;
  % set grid angle
  roms_angle=zeros(size(depth));

  % set masking
  mask_rho=ones(size(depth));
  dum=find(depth<-1.9);
  mask_rho(dum)=0;
  
  % set coriolis f
  f=zeros(size(depth))+4.988e-5; %20N
  
  % compute wave number
  
  L0=(9.81*TA.^2)./(2*pi);
  L=zeros(size(L0));
  
  for nn=1:100
      L=L0.*tanh((depth.*2*pi)./L);
  end
  
  k=(2*pi)./L;
  
  % compute wave celerity
  
  C=L./TA;
  
  % compute wave group celerity
  
  Cg=(C.*0.5).*(1+((k.*2).*depth)./sinh((k.*2).*depth));
  
  % compute refraction
  
  [mm,nn]=size(Cg);
  ang=zeros(size(Cg));
  
  for yy=1:mm
      ang(yy,:)=asin(C(yy,:)./C(yy,1)*sin((theta-270)*pi/180));
  end
  
  %ang=270+ang.*180/pi;
  
  % compute refraction coefficient
  Kr=zeros(size(Cg));
  
  for yy=1:mm
      Kr(yy,:)=(cos(ang(yy,:))/cos(ang(yy,1))).^0.5;
  end
 
  % compute shoaling
  
  Ks=zeros(size(Cg));
  
  for yy=1:mm
    Ks(yy,:)=(Cg(yy,1)./Cg(yy,:)).^0.5;
  end
  
  pcolor(Ks)
  shading flat
  colorbar 'vertical'
  
  ang=270+ang.*180/pi;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 INITIAL CONDITION DEFINITION                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_ini || make_InWave_bnd )  

  Nbins= 11;                               % number of directional bins considered in the simulation
  Bindirs_c = [260:10/(Nbins-1):270];   % center angles of the directional bins.
  Bindirs = [260-(10/(Nbins-1))/2:10/(Nbins-1):270+(10/(Nbins-1))/2];
  pd=ones(size(Bindirs)).*180./(Nbins-1);
  
end  

if (make_InWave_ini)  
    
  ini_file=strcat(filepath,'InWave_ini.nc');  % name of the initial file

  Ac=ones(Lm,Mm,Nbins).*0;
  Cx=ones(Lm-1,Mm,Nbins).*0;
  Cy=ones(Lm,Mm-1,Nbins).*0;
  Ct=ones(Lm,Mm,Nbins+1).*0;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 BOUNDARY CONDITION DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_bnd)

  bnd_file=strcat(filepath,'InWave_bnd.nc');  % name of the boundary file

  % Duration of the simulation and time increment for the boundaries
  
  dt= 1;                % time increment for the boundaries (seconds)
  drtn= 3600;         % total duration of the simulation (seconds)

  time=[0:dt:drtn];
  
  % Specify by 1 the open boundary: N E S W
  obc=[1 0 1 1];

  Nbins_bnd= 11;         % number of directional bins with energy at the boundaries
  dir_bnd= Bindirs_c;      % center angle (degrees of the bin containing energy at the boudaries
  dumd=find(dir_bnd==theta);
 
 % Specify number of directions at the boundaries (we have to specify at
 % least 2)
 
%   Nbins_bnd= 1;         % number of directional bins with energy at the boundaries
%   dir_bnd= 261;         % center angle (degrees of the bin containing energy at the boudaries)

if sum(ismember(dir_bnd,Bindirs_c)) ~= Nbins_bnd; 
  bin_error=1;
else
  bin_error=0;
end      

% SPECIFY CHARCATERISTICS OF BICHROMATIC WAVE GROUPS

a1=1.00;
a2=0.5;
f1=0.1;
T1=1./f1;
f2=0.12;  %0.115
T2=1./f2;
Tm=(T1+T2)./2;
TA=Tm;

close all

 % if obc(4)==1
      Ac_west=zeros(Mm,Nbins_bnd,length(time));
      for i=1:Mm
          time_r=time-((i-1)*20).*(tan(-(theta-270)*pi/180))./Cg(i,1);
          eta=a1*cos(2*pi*f1*time_r)+a2*cos(2*pi*f2*time_r);
          E=1/8*1025*9.81*(2.*abs(hilbert(eta))).^2;
          Ac1=E./(2*pi/Tm);
          Ac_west(i,dumd,:)=Ac1(:);
      end
 % end
  
 % if obc(1)==1
    Ac_north=zeros(Lm,Nbins_bnd,length(time));
    suma=0;
      for i=1:Lm
          suma=suma+((i-1)*20).*(tan(-(theta-270)*pi/180))/Cg(Mm,i);
          time_r=time-(Mm-1)*20/Cg(Mm,1)-suma;
          eta=Ks(Mm,i).*(a1*cos(2*pi*f1*time_r)+a2*cos(2*pi*f2*time_r));
          HH=2.*abs(hilbert(eta));
          HH(HH>0.78*depth(Mm,i))=0.78*depth(Mm,i);
          E=1/8*1025*9.81*(HH).^2;
          Ac1=E./(2*pi/Tm); 
          if (ang(Mm,i)==ang(Mm,i))
          dist=abs(dir_bnd-ang(Mm,i));
          dumd=find(dist==min(dist));
          end
          Ac1(isnan(Ac1)==1)=0.0;
          Ac_north(i,dumd,:)=Ac1(:);
      end
 % end

 % if obc(3)==1
      Ac_south=zeros(Lm,Nbins_bnd,length(time));
      suma=0;
      for i=1:Lm
          suma=suma+((i-1)*20).*(tan(-(theta-270)*pi/180))/Cg(1,i);
          time_r=time-(1-1)*20/Cg(1,1)-suma;
          eta=Ks(1,i).*(a1*cos(2*pi*f1*time_r)+a2*cos(2*pi*f2*time_r));
          HH=2.*abs(hilbert(eta));
          HH(HH>0.78*depth(1,i))=0.78*depth(1,i);
          E=1/8*1025*9.81*(2.*abs(hilbert(eta))).^2;
          Ac1=E./(2*pi/Tm); 
          if (ang(1,i)==ang(1,i))
          dist=abs(dir_bnd-ang(1,i));
          dumd=find(dist==min(dist));
          end
          Ac1(isnan(Ac1)==1)=0.0;
          Ac_south(i,dumd,:)=Ac1(:);
      end
%  end
    
%  if obc(2)==1
    Ac_east=zeros(length(time),Nbins_bnd,Mm);
%  end


  if obc(4)==1
    TA_west=TA;
  end
  if obc(2)==1
    TA_east=TA;
  end
  if obc(3)==1
    TA_south=TA;
  end
  if obc(1)==1
    TA_north=TA;
  end

end



