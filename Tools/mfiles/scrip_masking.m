% script scrip_masking.m
%
% adjusts the masking in the wrf scrip grid.
%
% jcw Sept 7, 2008
%     Feb 20, 2009
%

%%%%%%%%%% BEGIN user input   %%%%%%%%%%%%
%
%1) Enter name of roms netcdf grid file.
grid_file = 'joe_tc_coarse_grd.nc';

%2) Enter name of netcdf output file to use by scrip.
out_file = 'joe_tc_coarse_roms_scrip.nc';

%%%%%%%%%% END of user input   %%%%%%%%%%%%




%4) Set masking in roms2wrf file
% determine masking for wrf grid
% what cells are not in the roms grid?
ncclear

%get roms grid
ncload Z:\jwarner\Isabel\roms\forcings\USeast_grd3.nc
cornersx=[lon_rho(1,1) lon_rho(1,end) lon_rho(end,end) lon_rho(end,1)];
cornersy=[lat_rho(1,1) lat_rho(1,end) lat_rho(end,end) lat_rho(end,1)];

%get wrf grid
ncload z:\jwarner\Isabel\roms\run106\xlongu.nc
ncload z:\jwarner\Isabel\roms\run106\xlatv.nc
LP=size(XLONG_U,2);
MP=size(XLAT_V,1);
lon_rho_wrf=repmat(XLONG_U(1,:),MP,1);
lat_rho_wrf=repmat(XLAT_V(:,1),1,LP);


%find wrf points inside roms grid
wrfinroms=inpolygon(lon_rho_wrf,lat_rho_wrf,cornersx,cornersy);
%save wrfmask.mat

%cd Z:\jwarner\Isabel\roms\forcings
%now put this into the mask array
%load D:\data\Carolinas\hurricanes\Isabel\roms\run101\wrfmask.mat
ncload Z:\jwarner\Isabel\roms\forcings\ROMS_USeast_grd3_to_WRF_grd3_weights.nc
nc=netcdf('Z:\jwarner\Isabel\roms\forcings\ROMS_USeast_grd3_to_WRF_grd3_weights.nc','w')
count=0;
ztemp=zeros(LP*MP,1);
for jj=1:MP
  for ii=1:LP
    count=count+1;
    ztemp(count)=wrfinroms(jj,ii);
  end
end
nc{'dst_grid_imask'}(:) = ztemp(:);
clear ztemp
ncclose('nc')
