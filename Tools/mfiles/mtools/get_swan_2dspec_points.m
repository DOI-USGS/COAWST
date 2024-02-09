% get_swan_2dspec_points
%
% uses par and chd grids to determine 
% 2D spec points.
%
% jcw Apr 4, 2023
%

%0) cd to a working dir
cd E:\data3\Projects\NOPP_hurricanes_modeling\grids

%1) enter name of par grid
par_grid='GOMSAB_2km_ext_smooth.nc';

%2) enter name of chd grid
chd_grid='FL_Bbend2_CMX_h10min.nc';

%3) number of points to skip along each bndry
spec_res=25;

%4) tag for boundary points
bctag='FLBb';

%5)  time to start and duration, make this a string
out_time='20230827.000000 15 MIN';

%%%%%   end of user input   %%%%%%%%%%%%

%read in parent lon lat and mask
par_lon=ncread(par_grid,'lon_rho');
par_lat=ncread(par_grid,'lat_rho');
par_mask=ncread(par_grid,'mask_rho');

% read child grid and get boundary perimeter points
chd_lon=ncread(chd_grid,'lon_rho');
chd_lat=ncread(chd_grid,'lat_rho');
chd_mask=ncread(chd_grid,'mask_rho');

%now modify chd mask to block out any points that have parent mask=0
F = scatteredInterpolant(par_lon(:),par_lat(:),par_mask(:));
F.Method = 'nearest';
chd_mask2=F(chd_lon,chd_lat);
chd_mask2=chd_mask.*chd_mask2;
figure
hold on
pcolorjw(chd_lon,chd_lat,chd_mask2)
title('child mask from parent')

sr2=floor(spec_res/2);
%get southern points
ic=0;
for mm=2:spec_res:size(chd_lon,1)
  if(chd_mask2(mm,1)>0)
    ic=ic+1;
    bc_lon(ic)=chd_lon(mm,1);
    bc_lat(ic)=chd_lat(mm,1);
    xs(ic)=max(mm-sr2,1);
    xe(ic)=min(mm+sr2,size(chd_lon,1)-1);
    ys(ic)=1;
    ye(ic)=1;
  end
end
%get eastern points
for mm=2:spec_res:size(chd_lon,2)
  if(chd_mask2(end,mm)>0)
    ic=ic+1;
    bc_lon(ic)=chd_lon(end,mm);
    bc_lat(ic)=chd_lat(end,mm);
    xs(ic)=size(chd_lon,1)-1;
    xe(ic)=size(chd_lon,1)-1;
    ys(ic)=max(mm-sr2,1);
    ye(ic)=min(mm+sr2,size(chd_lon,2)-1);
  end
end
%get northern points
for mm=2:spec_res:size(chd_lon,1)
  if(chd_mask2(mm,end)>0)
    ic=ic+1;
    bc_lon(ic)=chd_lon(mm,end);
    bc_lat(ic)=chd_lat(mm,end);
    xs(ic)=max(mm-sr2,1);
    xe(ic)=min(mm+sr2,size(chd_lon,1)-1);
    ys(ic)=size(chd_lon,2)-1;
    ye(ic)=size(chd_lon,2)-1;
  end
end
%get western points
for mm=2:spec_res:size(chd_lon,2)
  if(chd_mask2(1,mm)>0)
    ic=ic+1;
    bc_lon(ic)=chd_lon(1,mm);
    bc_lat(ic)=chd_lat(1,mm);
    xs(ic)=1;
    xe(ic)=1;
    ys(ic)=max(mm-sr2,1);
    ye(ic)=min(mm+sr2,size(chd_lon,2)-1);
  end
end

%plot the points
plot(bc_lon,bc_lat,'r+')
for mm=1:ic
  text(bc_lon(mm),bc_lat(mm),num2str(mm))
end


%create text to add to parent SWAN input files
par_file='SWAN_par.txt';
fid=fopen(par_file,'w');
fprintf(fid,'& **   COPY THESE LINES TO SWAN PARENT INPUT FILE   ********\n');
%POINT 'MB1'   -86.2148   28.7289
%SPEC 'MB1' SPEC2D ABS 'MB1_2Dspec.txt' OUT 20181005.000000 1 HR
for mm=1:ic
  toadd=num2str(mm);
  nline=['POINT ''',bctag,toadd,''' ',num2str(bc_lon(mm)),'  ',num2str(bc_lat(mm))];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
%
  nline=['SPEC ''',bctag,toadd,'''  SPEC2D ABS  ', '''',bctag,toadd,'','spec.txt',''' OUT  ',out_time];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
end
fprintf(fid,'\n');
fclose(fid);

%create text to add to child SWAN input files
chd_file='SWAN_chd.txt';
fid=fopen(chd_file,'w');
fprintf(fid,'& **   COPY THESE LINES TO SWAN CHILD INPUT FILE   ********\n');
%BOUNDSPEC SEGMENT IJ     1    1     1    100 VARIABLE FILE 0 '../michael31/MB1_2Dspec.txt'
for mm=1:ic
  toadd=num2str(mm);
  locs=[num2str(xs(mm)),' ',num2str(ys(mm)),' ',num2str(xe(mm)),' ',num2str(ye(mm))];
  fname=['''',bctag,toadd,'spec.txt','''']
  nline=['BOUNDSPEC SEGMENT IJ ',locs,'  VARIABLE FILE 0  ', fname ];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
end
fprintf(fid,'\n');
fclose(fid);

disp('***    READ THE FILES  SWAN_chd.txt and SWAN_par.txt      ***')


