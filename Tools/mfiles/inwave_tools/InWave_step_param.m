% THIS FILE CONTAINS THE DEFINITIONS OF THE PARAMETERS NEEDED TO CREATE THE 
% INPUT FILES FOR THE INWAVE MODEL:
%
% InWave_grd.nc
% InWave_ini.nc
% InWave_bnd.nc
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                GRID AND BATHYMETRY DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid characteristics

Lm= 50;                % number of interior rho points in the xi direction
Mm= 20;                % number of interior rho points in the eta direction
dx=100;                % grid cell size in xi direction
dy=100;                % grid cell size in eta direction 

% Bathymetry characteristics

depth0= 15;            % water depth in the study domain (m)

% Wave characteristics

Nbins= 10;            % number of directional bins considered in the simulation
Bindirs = [0:360/Nbins:360-360/Nbins]; %center angles of the directional bins.

TA= 10;                % representative absolute wave period (sec)

% Specify by 1 the open boundaries

%    S E N W
obc=[0 0 0 1];

Nbins_bnd= 1;         % number of directional bins with energy at the boundaries
dir_bnd= 270;          % center angle (degrees of the bin containing energy at the boudaries)
idiot check = on.
  if sum(ismember(dir_bnd,Bindirs)) ~= Nbins_bnd; disp('Hey - you need to modify the bin dirs. One of the bnd dirs is not a bin dir."
disp(Bindirs, bnddslkdsalfj)

% Duration of the simulation and time increment for the boundaries

dt= 10;                % time increment for the boundaries (seconds)
drtn= 3600*3;          % total duration of the simulation (seconds)


% NC file names

grd_file='InWave_grd.nc';  % name of the grid file
ini_file='InWave_ini.nc';  % name of the initial file
bnd_file='InWave_bnd.nc';  % name of the boundary file

***********************************************************************************
dont go down here !!!  stay out !!!

% enter x and y coordinates of rho points
    x=[-dx/2:dx:dx*(Lm-1)];
    y=[-dy/2:dy:dy*(Mm-1)];

    x=repmat(x,length(y),1);
    y=repmat(y',1,length(x));
    
% set depth 
    depth=zeros(size(x))+depth0;

% set grid angle
    roms_angle=zeros(size(depth));

% set masking
    mask_rho=ones(size(depth));

% set coriolis f
    f=zeros(size(depth))+4.988e-5; %20N
dmde
dndx
ncload grid.nc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 INITIAL CONDITION DEFINITION                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=Lm-1;
M=Mm-1;

Ac=zeros(1,Nangle,Mm,Lm);
Cx=zeros(1,Nangle,Mm,L);
Cy=zeros(1,Nangle,M,Lm);
Ct=zeros(1,Nangle,Mm,Lm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 BOUNDARY CONDITION DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time=[0:Dt:Drtn];  

if obc(1)==1
Ac_north=zeros(length(time),Nangle_bnd,Lm);
end
if obc(2)==1
  Ac_east=zeros(length(time),Nangle_bnd,Mm);
end
if obc(3)==1
Ac_south=zeros(length(time),Nangle_bnd,Lm);
end
if obc(4)==1
Ac_west=zeros(length(time),Nangle_bnd,Mm);
Ac_west(60:60*3,1,:)=100;
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



