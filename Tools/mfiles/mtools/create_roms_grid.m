% create_roms_grid.m
%
% input:         See below.
% output:        Netcdf roms grid and ascii swan grids.
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
JOE_TC=0;
JOE_TC_coarse=0;
inlet_test=0;
inlet_test_diff=1;
MY_APP=0;

if (JOE_TC)
  %2) enter x and y coordinates of rho points
    ncellsx=199;  dx=12000;
    ncellsy=149;  dy=12000;
    xx=[-dx/2:dx:dx*(ncellsx-1)];
    yy=[-dy/2:dy:dy*(ncellsy-1)];
%
    x=repmat(xx,length(yy),1);
    y=repmat(yy',1,length(xx));
  %3) set depth 
    depth=zeros(size(x))+10;
    depth(:,50:69)=repmat(10+10*([50:69]-49),size(x,1),1);
    depth(:,70:89)=repmat(200+40*([70:89]-69),size(x,1),1);
    depth(:,90:end)=1000;
  %4) set grid angle
    roms_angle=zeros(size(depth));
  %5) set masking
    mask_rho=ones(size(depth));
    mask_rho(:,1:49)=0;
    mask_rho(:,end-3:end)=0;
    mask_rho(1:15,:)=0;
    mask_rho(end-3:end,:)=0;
  %6) set coriolis f
    f=zeros(size(depth))+4.988e-5; %20N
  %7) enter output file name
    fname='joe_tc_grd.nc';
elseif (JOE_TC_coarse)
  %2) enter x and y coordinates of rho points
    ncellsx=100;  dx=24000;
    ncellsy=75;   dy=24000;
    xx=[-dx/2:dx:dx*(ncellsx-1)]-48*111000;
    yy=[-dy/2:dy:dy*(ncellsy-1)]+12.4*111000;
  % 
    x=repmat(xx,length(yy),1);
    y=repmat(yy',1,length(xx));
  %3) set depth 
    depth=zeros(size(x))+10;
    depth(:,26:35)=repmat(10+20*([26:35]-25),size(x,1),1);
    depth(:,36:45)=repmat(200+80*([36:45]-35),size(x,1),1);
    depth(:,46:end)=1000;
  %4) set grid angle
    roms_angle=zeros(size(depth));
  %5) set masking
    mask_rho=ones(size(depth));
    mask_rho(:,1:25)=0;
    mask_rho(:,end-2:end)=0;
    mask_rho(1:8,:)=0;
    mask_rho(end-2:end,:)=0;
  %6) set coriolis f
    f=zeros(size(depth))+4.988e-5; %20N
  %7) enter output file name
    fname='joe_tc_coarse_grd.nc';
elseif (inlet_test)
  %2) enter x and y coordinates of rho points
    ncellsx=77;  dx=200;
    ncellsy=72;  dy=200;
    xx=[-100:dx:dx*(ncellsx-1)-100]; %16100];
    yy=[-100:dy:dy*(ncellsy-1)-100]; %15100];
  % 
    x=repmat(xx,length(yy),1);
    y=repmat(yy',1,length(xx));
  %3) set depth 
    depth=zeros(size(x))+4;
    zz=[4+0.0016*dy.*([37:72]-36)]';
    depth(37:72,:)=repmat(zz,1,length(xx));
  %4) set grid angle
    roms_angle=zeros(size(depth));
  %5) set masking
    mask_rho=ones(size(depth));
    mask_rho(1,:)=0;
    mask_rho(1:36,1)=0;
    mask_rho(1:36,end)=0;
    mask_rho(36,1:33)=0;
    mask_rho(36,45:end)=0;
  %6) set coriolis f
    f=zeros(size(depth))+4.988e-5;
  %7) enter output file name
    fname='inlet_test_grid.nc';
elseif (inlet_test_diff)
  %2) enter x and y coordinates of rho points
    ncellsx=87;  dx=200;
    ncellsy=82;  dy=200;
    xx=[-1100:dx:dx*(ncellsx-1)-1100]; %16100];
    yy=[-1100:dy:dy*(ncellsy-1)-1100]; %15100];
  % 
    x=repmat(xx,length(yy),1);
    y=repmat(yy',1,length(xx));
  %3) set depth 
    depth=zeros(size(x))+4;
    zz=[4+0.0016*dy.*([42:82]-41)]';
    depth(42:82,:)=repmat(zz,1,length(xx));
  %4) set grid angle
    roms_angle=zeros(size(depth));
  %5) set masking
    mask_rho=ones(size(depth));
    mask_rho(1:6,:)=0;
    mask_rho(1:40,1:6)=0;
    mask_rho(1:40,end-5:end)=0;
    mask_rho(41,1:38)=0;
    mask_rho(41,50:end)=0;
  %6) set coriolis f
    f=zeros(size(depth))+4.988e-5;
  %7) enter output file name
    fname='inlet_test_roms_bigger_grid.nc';
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

%to create swan grid then use:
  disp(['Created swan grid + bathy files:'])
  roms2swan(rho.x,rho.y,rho.depth,rho.mask);


