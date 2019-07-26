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
% updated to matlab netcdf 27Mar2012
%

%%%%%%%%%%%%% START OF USER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Select your app and then fill in steps 2-5 for that application.
% Required values are :
% x, y, depth, mask, and file name.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1) set your case here to = 1.
WET_DRY_DAM_BREAK=1;
MY_APP=0;

if (WET_DRY_DAM_BREAK)
  %2) enter x and y coordinates of rho points
    ncellsx=102;  dx=0.04;
    ncellsy=27;   dy=0.08;
    x=[-dx/2:dx:dx*(ncellsx-1)];
    y=[-dy/2:dy:dy*(ncellsy-1)];
  % 
    x=repmat(x,length(y),1)';
    y=repmat(y',1,length(x))';
  %3) set depth 
    depth=x*0+0.001;
  %4) set masking
    mask_rho=ones(size(depth));
    mask_rho(1,:)=0;
    mask_rho(1:26,1)=0;
    mask_rho(1:26,end)=0;
    mask_rho(26,1:11)=0;
    mask_rho(26,17:end)=0;
  %5) enter output file name
    fname='wetdry_dam_break_grd.nc';

elseif (MY_APP)
  disp('set MY_APP=1 and then put your stuff in here')
end

%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

%create roms grid
  rho.x=x;
  rho.y=y;  
  rho.depth=depth;
  rho.mask = mask_rho
  projection='mercator';
  spherical='F';
  save temp_jcw33.mat
  eval(['mat2roms_mw(''temp_jcw33.mat'',''',fname,''');'])
  !del temp_jcw33.mat
  disp(['Created roms grid -->   ',fname])


