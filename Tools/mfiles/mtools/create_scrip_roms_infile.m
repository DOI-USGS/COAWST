function create_scrip_roms_infile(roms_grid_file,mask,out_file);

%jcw 30July2014 from srcip_roms to make it a function.

%load the roms grid
%netcdf_load(roms_grid_file)
h=ncread(roms_grid_file,'h');
x_psi=ncread(roms_grid_file,'lon_psi');
y_psi=ncread(roms_grid_file,'lat_psi');
%mask_rho=ncread(roms_grid_file,'mask_rho');
mask_rho=mask;
spherical=ncread(roms_grid_file,'spherical');

[LP,MP]=size(h);
%gridsize=LP*MP;

if ((spherical=='F')||(spherical=='f'))||(spherical==0)
%  lon_rho=x_rho;
%  lat_rho=y_rho;
%  lon_psi=x_psi;
%  lat_psi=y_psi;
   x_psi=ncread(roms_grid_file,'x_psi');
   y_psi=ncread(roms_grid_file,'y_psi');
   x_rho=ncread(roms_grid_file,'x_rho');
   y_rho=ncread(roms_grid_file,'y_rho');
   projection='mercator';
   m_proj(projection);
   [lon_rho, lat_rho] = m_xy2ll(x_rho/6371000, y_rho/6371000);   % Degrees.
   [lon_psi, lat_psi] = m_xy2ll(x_psi/6371000, y_psi/6371000);   % Degrees.
else
  lon_psi=ncread(roms_grid_file,'lon_psi');
  lat_psi=ncread(roms_grid_file,'lat_psi');
  lon_rho=ncread(roms_grid_file,'lon_rho');
  lat_rho=ncread(roms_grid_file,'lat_rho');
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

