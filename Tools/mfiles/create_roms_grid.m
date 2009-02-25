% create_roms_grid.m
%
% input:         x = array x values
%                y = array of y values
%                depth = array of depth values
%                angle = assumed to = 0.
%                ncfile_name = name for the grid to be produced.
% output:        netcdf roms grid.
%
% jcwarner March 3, 2007
%

%%%%%%%%%%%%% START OF USER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% select your app and then fill in steps 2-7 for that application %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) set your case here to = 1.
WET_DRY_SLOPE_CHAN=0;
JOE_TC=1;
MY_APP=0;

if (WET_DRY_SLOPE_CHAN)
  %2) enter x and y coordinates of rho points
    ncellsx=101;  dx=250;
    ncellsy=6;    dy=200;
    x=[-125:dx:25125];
    y=[-100:dy:1100];
  % 
    x=repmat(x,length(y),1);
    y=repmat(y',1,length(x));
  %3) set depth 
    depth=10*x/25125;
  %4) set grid angle
    roms_angle=zeros(size(depth));
  %5) set masking
    mask_rho=ones(size(depth));
  %6) set coriolis f
    f=zeros(size(depth));
  %7) enter output file name
    fname='ocean_wetdry_slope_chan_grd.nc';
elseif (JOE_TC)
  %2) enter x and y coordinates of rho points
    ncellsx=200;  dx=12000;
    ncellsy=150;  dy=12000;
    x=[-dx/2:dx:dx*(ncellsx-1)];
    y=[-dy/2:dy:dy*(ncellsy-1)];
  % 
    x=repmat(x,length(y),1);
    y=repmat(y',1,length(x));
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


