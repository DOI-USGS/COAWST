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

% 6) Enter name of WW3 unstructured msh file to be created:
ww3_msh_file='ww3_sandy_grid.msh';

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
%h(h<(min_depth))=9999;
for index1 = 1:MP;
    for index2 = 1:LP;
       fprintf(fid,'   ');
       fprintf(fid,'%12.2f',h(index2,index1));
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
ww3_mask(:,1)=2;
ww3_mask(:,end)=2;
ww3_mask(1,:)=2;
ww3_mask(end,:)=2;
%
% now remove land and set these points as 0
%
for index1 = 2:MP-1;
    for index2 = 2:LP-1;
        if (mask_rho(index2,index1)==0)
          ww3_mask(index2,index1)=0;
        end
        if (h(index2,index1)<min_depth)
%          ww3_mask(index2,index1)=0;
        end
    end
end
%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  here we create the unstructured mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
  lat_rho=lat_rho+0.00001*rand(1);
  lon_rho=lon_rho+0.00001*rand(1);
end

%
% set grid increments
%
xs=1;xe=size(h,1);
ys=1;ye=size(h,2);
h=h(xs:xe,ys:ye);
h(mask_rho(xs:xe,ys:ye)==0)=0;
lon_rho=lon_rho(xs:xe,ys:ye);
lat_rho=lat_rho(xs:xe,ys:ye);
[LP, MP]=size(h);

%open output file for writing
fid=fopen(ww3_msh_file,'w');

%write the header stuff
fprintf(fid,'%11s\n','$MeshFormat');
fprintf(fid,'%5s\n','2 0 8');
fprintf(fid,'%14s\n','$EndMeshFormat');
fprintf(fid,'%6s\n','$Nodes');

%compute number of nodes
if (spherical=='T' || spherical=='t' || spherical==1)
  x_full_grid=lon_rho;
  y_full_grid=lat_rho;
else
  x_full_grid=x_rho;
  y_full_grid=y_rho;
end
%
%  "The Nodes section starts with the total number of nodes (grid points),
%  and then lists the details of each one, namely its index, 
%  longitude (x), latitude (y) and elevation: "
%
%  write num nodes
Numnodes=size(x_full_grid,1)*size(x_full_grid,2);
fprintf(fid,'%d\n',Numnodes);

%now write out each node: number lon lat depth
nodesx=x_full_grid(:);
nodesy=y_full_grid(:);
nodesh=h(:);
for mm=1:Numnodes
  fprintf(fid,'%d  %f  %f  %f\n',[mm nodesx(mm) nodesy(mm) h(mm)]);
end
fprintf(fid,'%6s\n','$EndNodes');

%compute elements
fprintf(fid,'%6s\n','$Elements');

%  write num elems
Numelems=(size(x_full_grid,1)-1)*(size(x_full_grid,2)-1)*2+size(x_full_grid,1)*2+(size(x_full_grid,2)-2)*2;
fprintf(fid,'%d\n',Numelems);

% "Finally, the Elements section provides information on the connectivity
% of these nodes, describing the mesh elements and boundary locations. The 
% first number in this section is the total number of elements, followed by 
% a listing of all the elements and their characteristics. These include: 
% element index, element type (2 = 3-node triangle; 15 = 1-node point, 
% used to describe boundary points), total number of tags describing this element,
% the list of these tags (e.g. used to specify boundary conditions along a given grid edge), 
% and the list of nodal indices forming the corner points of the element 
% (=3 for triangular elements; =1 for boundary points). "

%lets write out the boundaries first, all 4 sides
%
%south
count=0;
for mm=1:size(x_full_grid,1)
  count=count+1;
  fprintf(fid,'%d  %d  %d  %d  %d  %d\n',[count 15  2  0  0  count]);
end
%
%north
X=size(x_full_grid,1)*(size(x_full_grid,2)-1);
for mm=1:size(x_full_grid,1)
  count=count+1;
  X=X+1;
  fprintf(fid,'%d  %d  %d  %d  %d  %d\n',[count 15  2  0  0  X]);
end
%
%west, but not SW or NW corners, already included.
for mm=2:size(x_full_grid,2)-1
  count=count+1;
  X=size(x_full_grid,1)*(mm-1)+1;
  fprintf(fid,'%d  %d  %d  %d  %d  %d\n',[count 15  2  0  0  X]);
end
%
%east, but not NE or SE corners, already included.
for mm=2:size(x_full_grid,2)-1
  count=count+1;
  X=size(x_full_grid,1)*mm;
  fprintf(fid,'%d  %d  %d  %d  %d  %d\n',[count 15  2  0  0  X]);
end
%
% now the working 3 point elements. start LL.  CW.
%
count2=0;
for jj=1:size(x_full_grid,2)-1
  for ii=1:size(x_full_grid,1)-1
% lower left
    count=count+1;
    count2=count2+1;
    X1=size(x_full_grid,1)*(jj-1)+ii;
    X2=size(x_full_grid,1)*jj+ii;
    X3=size(x_full_grid,1)*(jj-1)+ii+1;
    fprintf(fid,'%d  %d  %d  %d  %d  %d  %d  %d  %d\n',[count 2  3  0  count2  0  X1  X2  X3]);
% upper right
    count=count+1;
    count2=count2+1;
    X1=size(x_full_grid,1)*(jj-1)+ii+1;
    X2=size(x_full_grid,1)*jj+ii;
    X3=X2+1;
    fprintf(fid,'%d  %d  %d  %d  %d  %d  %d  %d  %d\n',[count 2  3  0  count2  0  X1  X2  X3]);
  end
end
%
% close
%
fprintf(fid,'%6s\n','$EndElements');
fclose(fid)



%%%%%%%%%%%
% example here
%
% $MeshFormat
% 2 0 8
% $EndMeshFormat
% $Nodes
% 8639
% 1  -87.4999338083   30.2829507916     0.0000000000
% 2  -87.4999324079   30.2668530854     7.4275088310
% 3  -87.4999308308   30.2481519873     9.3883380890
% 4  -87.4999291947   30.2280422425    12.7037878036
% ...
% 8636  -90.4291892729   30.1571185245     0.0000000000
% 8637  -90.4285705858   30.1436015971     0.4048240483
% 8638  -89.4482000000   30.0026000000     0.0676150768
% 8639  -89.7497900000   30.1078600000     0.0476770701
% $EndNodes
% $Elements
% 16245
% 1  15  2  0  0  1
% 2  15  2  0  0  2
% 3  15  2  0  0  3
% 4  15  2  0  0  4
% ...
% 
% 126  15  2  0  0  126
% 127  15  2  0  0  127
% 128  2  3  0  1  0  1  128  2
% 129  2  3  0  2  0  129  2  128
% 130  2  3  0  3  0  3  2  129
% 131  2  3  0  4  0  130  3  129
% 132  2  3  0  5  0  3  130  131
% 133  2  3  0  6  0  131  4  3
% 134  2  3  0  7  0  131  132  4
% 135  2  3  0  8  0  4  132  5
% 
% 16243  2  3  0  16116  0  8638  7330  7299
% 16244  2  3  0  16117  0  7909  7934  8639
% 16245  2  3  0  16118  0  7886  7909  8639
% $EndElements
% 
% 


