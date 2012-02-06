% create_roms_xygrid.m
%
% input:         Grid locations (see below).
% output:        Netcdf roms grid and ascii swan grids.
%
% This is intended to be a simple roms grid, typically
% rectilinear, in x-y space, with spherical=F. 
% For more complicated grids use seagrid or some other tool.
% To create a roms grid from wrf, use wrf2roms_mw.m
%
% Note: This uses the NATIVE MATLAB interface.
%
%  help netcdf.defVar provides:
% " Because MATLAB uses FORTRAN-style ordering, the fastest-varying 
%   dimension comes first and the slowest comes last. Any unlimited
%   dimension is therefore last in the list of dimension IDs. This 
%   ordering is the reverse of that found in the C API."
%
% The variable order is (x,y,z,t), different than previous methods
%
% jcwarner 29Jan2009: original
% jcwarner 31Jan2012: updated to use newer version of mat2roms_mw.m
%

%%%%%%%%%%%%% START OF USER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Select your app and then fill in steps 2-7 for that application.
% Required values are :
%
% - rho point locations/values of x, y, depth, and mask_rho.
% - provide output netcdf file name.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1) set your case here to = 1.
inlet_test=1;
inlet_test_diff=0;
MY_APP=0;

if (inlet_test)
  %2) enter x and y coordinates of rho points
    ncellsx=77;  dx=200;
    ncellsy=72;  dy=200;
    xx=[-100:dx:dx*(ncellsx-1)-100];
    yy=[-100:dy:dy*(ncellsy-1)-100];
  % 
    x=repmat(xx',1,length(yy));
    y=repmat(yy,length(xx),1);
  %3) set depth 
    depth=zeros(size(x))+4;
    zz=[4+0.0016*dy.*([37:72]-36)];
    depth(:,37:72)=repmat(zz,length(xx),1);
  %4) set masking
    mask_rho=ones(size(depth));
    mask_rho(:,1)=0;
    mask_rho(1,1:36)=0;
    mask_rho(end,1:36)=0;
    mask_rho(1:33,36)=0;
    mask_rho(45:end,36)=0;
  %5) enter output file name
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
  %4) set masking
    mask_rho=ones(size(depth));
    mask_rho(1:6,:)=0;
    mask_rho(1:40,1:6)=0;
    mask_rho(1:40,end-5:end)=0;
    mask_rho(41,1:38)=0;
    mask_rho(41,50:end)=0;
  %5) enter output file name
    fname='inlet_test_roms_bigger_grid.nc';

elseif (MY_APP)
  disp('set MY_APP=1 and then put your stuff in here')
end

% for all use:
    spherical='F';
    projection='mercator';

%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

%create roms grid
  rho.x=x;
  rho.y=y;  
  rho.depth=depth;
  rho.mask = mask_rho;
  save temp_jcw33.mat
  eval(['mat2roms_mw(''temp_jcw33.mat'',''',fname,''');'])
  !del temp_jcw33.mat
  disp(['Created roms grid -->   ',fname])

%to create swan grid then use:
  disp(['Created swan grid + bathy files:'])
  roms2swan(rho.x,rho.y,rho.depth,rho.mask);


