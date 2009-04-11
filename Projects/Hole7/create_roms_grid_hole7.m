% create_roms_grid.m
%
% input:         See below.
% output:        Netcdf roms grid.
%
% This is intended to be a simple roms grid, typically
% rectilinear. For more complicated grids use seagrid of
% something else.
%
% jcwarner 21Feb2009
%

%%%%%%%%%%%%% START OF USER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Select your app and then fill in steps 2-7 for that application.
% Required values are :
% x, y, dx, dy, depth, angle, mask, f, and file name.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1) set your case here to = 1.
HOLE7=1;

if (HOLE7)
  %2) enter x and y coordinates of rho points
    ncellsx=151;      dx=20;
%   ncellsy=196+2*40; dy=80;
    ncellsy=196; dy=80;
    xx=[-dx/2:dx:dx*(ncellsx-1)];
    yy=[-dy/2:dy:dy*(ncellsy-1)];
    yy=yy+3200;
%
    x=repmat(xx,length(yy),1);
    y=repmat(yy',1,length(xx));
  %3) set depth 
    depth=zeros(size(x))+10;
  %4) set grid angle
    roms_angle=zeros(size(depth));
  %5) set masking
    mask_rho=ones(size(depth));
  %6) set coriolis f
    f=zeros(size(depth))+4.988e-5; %20N
  %7) enter output file name
%   fname='hole7_wave_grid.nc';
    fname='hole7_ocean_grid.nc';
elseif (MY_APP)
  disp('set MY_APP=1 and then put your stuff in here')
end

%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

%create roms grid
  rho.x=x;
  rho.y=y;  
  rho.dx=dx;
  rho.dy=dy;
  rho.depth=depth;
  rho.angle = roms_angle;
  rho.mask = mask_rho
  rho.f = f;
  projection='mercator';
  save temp_jcw33.mat
  eval(['mat2roms_jcw(''temp_jcw33.mat'',''',fname,''');'])
  !del temp_jcw33.mat
  disp(['Created roms grid -->   ',fname])


%for ocean grid
plot_hole07_bathy
ncload hole7_ocean_grid.nc
h=depth_array_mod;
nc=netcdf('hole7_ocean_grid.nc','w');
nc{'h'}(:)=h(:);
ncclear
  

%for wave grid
ncload hole7_ocean_grid.nc
hocean=h;
ncload hole7_wave_grid.nc
h(1+40:196+40,:)=hocean;
for mm=1:40
  h(mm,:)=h(41,:);
end
for mm=237:276
  h(mm,:)=h(236,:);
end
nc=netcdf('hole7_wave_grid.nc','w');
nc{'h'}(:)=h(:);
ncclear
disp(['Created swan grid + bathy files:'])
roms2swan(x_rho,y_rho,h,mask_rho);


