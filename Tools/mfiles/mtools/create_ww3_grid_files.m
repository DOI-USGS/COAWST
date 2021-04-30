% create_ww3_grid_files
%
% This mfile will create WAVEWATCH III grid files of:
% xcoord, ycoord, bathy, and mapsta.
%
% This m file is set up to work on the Projects/Sandy test case. IT can be 
% modified for other applications.
%
% jcwarner 13Nov2017
%

% 1) Enter name of ROMS grid, or provide the varaibles of 
% lon_rho, lat_rho, h, and mask_rho:
roms_grid='Sandy_roms_grid.nc';
% or 
% lon_rho= , lat_rho= , mask_rho= , h= .

% 2) Enter name of WW3 xcoord file to be created:
ww3_xcoord_file='ww3_sandy_xcoord.dat';

% 3) Enter name of WW3 ycoord file to be created:
ww3_ycoord_file='ww3_sandy_ycoord.dat';

% 4) Enter name of WW3 bathy file to be created:
ww3_bath_file='ww3_sandy_bathy.bot';
min_depth=5.1;

% 5) Enter name of WW3 mask file to be created:
ww3_mask_file='ww3_sandy_mapsta.inp';

%%%%%%%%%%%%%%%%%%    END OF USER INPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) Start by getting roms grid. For this test we will create a WW3 grid that is the 
% same as ROMS.
%
netcdf_load(roms_grid)
[LP, MP]=size(h);
%
% 2) Create WW3 x coords
%
fid = fopen(ww3_xcoord_file,'w');
for index1 = 1:MP;
    for index2 = 1:LP;
       fprintf(fid,'   ');
       if (spherical=='T' || spherical=='t' || spherical==1)
         fprintf(fid,'%12.8f',lon_rho(index2,index1));
       else
         fprintf(fid,'%12.8f',x_rho(index2,index1));
       end
    end
    fprintf(fid,'\n');
end
fclose(fid);
%
% 3) Create WW3 y coords
%
fid = fopen(ww3_ycoord_file,'w');
for index1 = 1:MP;
    for index2 = 1:LP;
       fprintf(fid,'   ');
       if (spherical=='T' || spherical=='t' || spherical==1)
         fprintf(fid,'%12.8f',lat_rho(index2,index1));
       else
         fprintf(fid,'%12.8f',y_rho(index2,index1));
       end
    end
    fprintf(fid,'\n');
end
fclose(fid);
%
% 4) Create ww3 bathy
%
fid = fopen(ww3_bath_file,'w');
h(h<(min_depth))=9999;
for index1 = 1:MP;
    for index2 = 1:LP;
       fprintf(fid,'   ');
       fprintf(fid,'%12.8f',h(index2,index1));
    end
    fprintf(fid,'\n');
end
fclose(fid);
%
% 5) Create ww3 masking file.
%  Start with mask of all ones to set as all ocean
%$ The legend for the input map is :
%$    0 : Land point.
%$    1 : Regular sea point.
%$    2 : Active boundary point.
%$    3 : Point excluded from grid.
%
ww3_mask=ones(LP,MP);
%
% now make boudaries = 2
%
%ww3_mask(:,1)=2;
%ww3_mask(:,end)=2;
%ww3_mask(1,:)=2;
%ww3_mask(end,:)=2;
%
% now remove land and set these points as 0
%
zz=find(mask_rho==0);
ww3_mask(zz)=0;

% now write this out
%
fid = fopen(ww3_mask_file,'w');
for index1 = 1:MP;
    for index2 = 1:LP;
%      fprintf(fid,'   ');
       fprintf(fid,'%6i',ww3_mask(index2,index1));
    end
    fprintf(fid,'\n');
end
fclose(fid);


