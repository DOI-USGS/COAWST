%create_swan_refined_bc_points
%
%jcwarner Jan 20,, 2009
%
% determine the boudnary points for swan refined grids
%

%1) enter name of netcdf grid file
%ncfile='CH_shoals_grd4.nc';
ncfile='KH_CH_grd4_min5.nc';

%2) enter number of points / side
incx=12;
incy=12;

%%%%%%%%%%%  end of user input %%%%%%%%%%%%%%%%%%
eval(['ncload ',ncfile,' ']);

figure
pcolorjw(lon_rho,lat_rho,h);colorbar
hold on

%get bndry segments - lat lon style
%west
bndry_segs1=[lon_rho(1:incy:end-incy,1) lat_rho(1:incy:end-incy,1)];
bndry_segs1=[bndry_segs1; lon_rho(end,1) lat_rho(end,1)];
plot(bndry_segs1(:,1),bndry_segs1(:,2),'r+')
%east
bndry_segs2=[lon_rho(1:incy:end-incy,end) lat_rho(1:incy:end-incy,end)];
bndry_segs2=[bndry_segs2; lon_rho(end,end) lat_rho(end,end)];
plot(bndry_segs2(:,1),bndry_segs2(:,2),'k+')
%south
bndry_segs3=[lon_rho(1,1:incx:end-incx)' lat_rho(1,1:incx:end-incx)'];
bndry_segs3=[bndry_segs3; lon_rho(1,end) lat_rho(1,end)];
plot(bndry_segs3(:,1),bndry_segs3(:,2),'b+')
%north
bndry_segs4=[lon_rho(end,1:incx:end-incx)' lat_rho(end,1:incx:end-incx)'];
bndry_segs4=[bndry_segs4; lon_rho(end,end) lat_rho(end,end)];
plot(bndry_segs4(:,1),bndry_segs4(:,2),'m+')
%all 4
bndry_segs=[bndry_segs1 ; bndry_segs2 ; bndry_segs3 ; bndry_segs4];

%get bndry segments - IJ indices style
bndry_segsy=[0:incy:size(lon_rho,1)-incy size(lon_rho,1)-1]';
bndry_segsx=[0:incx:size(lon_rho,2)-incx size(lon_rho,2)-1]';
bndry_segsxy=[ [bndry_segsy(1:end-1)*0 bndry_segsy(1:end-1)   ] ...
               [bndry_segsy(2:end  )*0 bndry_segsy(2:end  )   ];
               [bndry_segsy(1:end-1)*0+bndry_segsx(end)  bndry_segsy(1:end-1)] ...
               [bndry_segsy(1:end-1)*0+bndry_segsx(end)  bndry_segsy(2:end  )];

               [bndry_segsx(1:end-1) bndry_segsx(1:end-1)*0   ] ...
               [bndry_segsx(2:end  ) bndry_segsx(2:end  )*0   ];

               [bndry_segsx(1:end-1) bndry_segsx(1:end-1)*0+bndry_segsy(end)   ] ...
               [bndry_segsx(2:end  ) bndry_segsx(2:end  )*0+bndry_segsy(end)   ];
              ];

% Do not write out the lat lons. SWAN may not interpret the points to be
% exactly on the bndry and give an error. So just use the IJ coords to be 
% more exact.

%& Boundary files  ****************************************
%BOUND SHAPESPEC JONSWAP 3.3 PEAK DSPR DEGREES
%'BOUNDSPEC SEGMENT IJ '  'CONSTANT PAR 0.1 10.0 0. 20.'
disp('cut the following lines and paste into INPUT for --- CHILD --- grid.')
disp('BOUND SHAPESPEC JONSWAP 3.3 PEAK DSPR DEGREES')
for mm=1:size(bndry_segsxy,1)
  disp(['BOUNDSPEC SEGMENT IJ ',num2str(bndry_segsxy(mm,:)),  '  CONSTANT PAR 0.1 10.0 0. 20.'])
end

%**************************************************************************
%now get bndry points
temp1=[bndry_segs1(2:end,:)-bndry_segs1(1:end-1,:)]/2+bndry_segs1(1:end-1,:);
bndry_pts1=temp1;
plot(bndry_pts1(:,1),bndry_pts1(:,2),'ro')

temp1=[bndry_segs2(2:end,:)-bndry_segs2(1:end-1,:)]/2+bndry_segs2(1:end-1,:);
bndry_pts2=temp1;
plot(bndry_pts2(:,1),bndry_pts2(:,2),'ko')

temp1=[bndry_segs3(2:end,:)-bndry_segs3(1:end-1,:)]/2+bndry_segs3(1:end-1,:);
bndry_pts3=temp1;
plot(bndry_pts3(:,1),bndry_pts3(:,2),'bo')

temp1=[bndry_segs4(2:end,:)-bndry_segs4(1:end-1,:)]/2+bndry_segs4(1:end-1,:);
bndry_pts4=temp1;
plot(bndry_pts4(:,1),bndry_pts4(:,2),'ro')

bndry_pts=[bndry_pts1 ; bndry_pts2 ; bndry_pts3 ; bndry_pts4];

disp('cut the following lines and paste into INPUT for --- PARENT --- grid.')
for mm=1:size(bndry_pts,1)
  disp(['POINTS ''point',num2str(mm), '''  ',num2str(bndry_pts(mm,:))])
end
for mm=1:size(bndry_pts,1)
  disp(['SPECOUT ''point',num2str(mm), ''' SPEC2D ''point',num2str(mm),'.spc2d'' OUTPUT 20000101.000000 10 MIN'])
end



