% script scrip_roms.m (formerly roms2scrip.m)
%
% convert roms grid to scrip input format
%
% jcw Sept 7, 2008
% jcw Feb 5, 2012 convert to matlab netcdf interface
%

clear
%%%%%%%%%% BEGIN user input   %%%%%%%%%%%%
%
%1) Enter name of roms netcdf grid file.
%grid_file = 'joe_tc_coarse_grd.nc';
grid_file = 'inlet_test_grid.nc';
%grid_file = 'inlet_test_grid_ref3.nc';

%2) Enter name of netcdf output file to use by scrip.
%out_file = 'joe_tc_coarse_roms_scrip.nc';
out_file = 'inlet_test_roms_scrip.nc';
%out_file = 'inlet_test_roms_ref3_scrip.nc';

%%%%%%%%%% END of user input   %%%%%%%%%%%%

%load the roms grid
netcdf_load(grid_file)

[LP,MP]=size(h);
%gridsize=LP*MP;

if ((spherical=='F')||(spherical=='f'))||(spherical==0)
%  lon_rho=x_rho;
%  lat_rho=y_rho;
%  lon_psi=x_psi;
%  lat_psi=y_psi;
   projection='mercator';
   m_proj(projection);
   [lon_rho, lat_rho] = m_xy2ll(x_rho/6371000, y_rho/6371000);   % Degrees.
   [lon_psi, lat_psi] = m_xy2ll(x_psi/6371000, y_psi/6371000);   % Degrees.
end

%create a full set of psi points
[x_full_grid,y_full_grid]=create_extra_rho_grid(lon_psi,lat_psi);

%create the srcip netcdf grid file
nc = netcdf.create(out_file, 'clobber');

%% Global attributes:
disp(' ## Defining Global Attributes...')

netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title','Scrip file');
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history',['Created by ', mfilename ', on ', datestr(now)]);

L = LP-1;       % The psi dimension.
M = MP-1;       % The psi dimension.

disp(' ## Defining Dimensions...')

grid_size = netcdf.defDim(nc,'grid_size',LP*MP);
grid_corners = netcdf.defDim(nc,'grid_corners',4);
grid_rank = netcdf.defDim(nc,'grid_rank',2);

disp(' ## Defining Variables')

v1 = netcdf.defVar(nc,'grid_dims','short',grid_rank);
netcdf.putAtt(nc,v1,'long_name','grid dimensions')
netcdf.putAtt(nc,v1,'units','---')

v2 = netcdf.defVar(nc,'grid_imask','short',grid_size);
netcdf.putAtt(nc,v2,'long_name','grid masking')
netcdf.putAtt(nc,v2,'units','---')

v3 = netcdf.defVar(nc,'grid_center_lat','double',grid_size);
netcdf.putAtt(nc,v3,'units','radians')

v4 = netcdf.defVar(nc,'grid_center_lon','double',grid_size);
netcdf.putAtt(nc,v4,'units','radians')

v5 = netcdf.defVar(nc,'grid_corner_lat','double',[grid_corners grid_size]);
netcdf.putAtt(nc,v5,'units','radians')

v6 = netcdf.defVar(nc,'grid_corner_lon','double',[grid_corners grid_size]);
netcdf.putAtt(nc,v6,'units','radians')

netcdf.endDef(nc)

%now fill that netcdf file
netcdf.putVar(nc,v1, [MP, LP]);

disp('step 0/4, filling grid mask')
%get grid mask

mask_rho=mask_rho*0;
Istr=25; Iend=54; Jstr=41; Jend=56;
mask_rho(Istr:Iend,Jstr:Jend)=1;
%mask_rho(:,1)=0;
%mask_rho(:,end)=0;
%mask_rho(1,:)=0;
%mask_rho(end,:)=0;

count=0;
for jj=1:MP
  for ii=1:LP
    count=count+1;
    ztemp(count)=mask_rho(ii,jj);
  end
end
netcdf.putVar(nc,v2,ztemp);
clear ztemp

scale=pi/180;

%get grid centers
disp('step 1/4, filling grid lat centers')
ztemp=reshape(lat_rho,MP*LP,1);
netcdf.putVar(nc,v3,ztemp*scale);
clear ztemp

disp('step 2/4, filling grid lon centers')
ztemp=reshape(lon_rho,MP*LP,1);
netcdf.putVar(nc,v4,ztemp*scale);
clear ztemp

%get grid corners, counterclockwise
disp('step 3/4, filling grid lat corners')
c1=reshape(y_full_grid(1:LP,1:MP),MP*LP,1);
c2=reshape(y_full_grid(2:LP+1,1:MP),MP*LP,1);
c3=reshape(y_full_grid(2:LP+1,2:MP+1),MP*LP,1);
c4=reshape(y_full_grid(1:LP,2:MP+1),MP*LP,1);
ztemp(:,:)=[c1 c2 c3 c4].';
netcdf.putVar(nc,v5,ztemp*scale);
clear ztemp

disp('step 4/4, filling grid lon corners')
c1=reshape(x_full_grid(1:LP,1:MP),MP*LP,1);
c2=reshape(x_full_grid(2:LP+1,1:MP),MP*LP,1);
c3=reshape(x_full_grid(2:LP+1,2:MP+1),MP*LP,1);
c4=reshape(x_full_grid(1:LP,2:MP+1),MP*LP,1);
ztemp(:,:)=[c1 c2 c3 c4].';
netcdf.putVar(nc,v6,ztemp*scale);
clear ztemp

%close the file.
netcdf.close(nc)

