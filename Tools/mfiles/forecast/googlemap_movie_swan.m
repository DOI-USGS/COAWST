%
% create movie of COAWST output
%
% BNA on 3/26/09

cd D:\Temporary\COAWST\SWAN\run\pngfiles

%enter name of output file
nc=netcdf('D:\Temporary\COAWST\SWAN\run\hsig.nc');
ncw=netcdf('D:\Temporary\COAWST\SWAN\run\wind.nc');
ncg=netcdf('D:\Temporary\COAWST\SWAN\forcing\USeast_grd3.nc');

frame=1;

%big do loop
ocean_time=nc{'wave_time'}(:);%/3600/24;
angle_r=ncg{'angle'}(:);
lon_rho=ncg{'lon_rho'}(:);
lat_rho=ncg{'lat_rho'}(:);
mask_rho=ncg{'mask_rho'}(:);
mask_rho(isnan(mask_rho))=1;
mask_rho(mask_rho==0)=nan;
minh=0;
maxh=0;
for ct=1:28;%length(ocean_time);
    minh_nw=min(min(nc{'Hwave'}(ct,:,:)));
    maxh_nw=max(max(nc{'Hwave'}(ct,:,:)));
    if minh>minh_nw
        minh=minh_nw;
    end
    if maxh<maxh_nw
        maxh=maxh_nw;
    end
end
close(ncg)
% figure;hold on;
% title_text=('Significant Wave Height');
% th=text(-0.05,0.05,title_text,'Rotation',90,'FontSize',20);
hc=colorbar('FontSize',20,'Location','West');
caxis([minh maxh+1]);
axis off
eval(['print -dpng -r0 D:/Temporary/COAWST/jet_colorbar.png;']);
%eval(['print -dpng -r0 Z:/testing/images/jet_colorbar.png;']);
crop('D:/Temporary/COAWST/jet_colorbar.png');
cd D:\Temporary\COAWST\
!scp -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" jet_colorbar.png barmstrong@capecodder.er.usgs.gov:/Volumes/web/user/cccp/public/images/
cd D:\Temporary\COAWST\SWAN\run\pngfiles

close all;
%[lon_rho_rot lat_rho_rot]=rot(lon_rho,lat_rho,angle_r);
figure
set(gcf,'position',[100 255 810 800])

for tidx=1:28%length(ocean_time)%length(ocean_time):length(ocean_time)

    Hwave=nc{'Hwave'}(tidx,:,:);
    Windx=ncw{'Windx'}(tidx,1:20:end,1:20:end);
    Windy=ncw{'Windy'}(tidx,1:20:end,1:20:end);
    clf
    hold on;
    pc_hwave=pcolorjw(lon_rho,lat_rho,Hwave.*mask_rho);
    pc_wind=quiver(lon_rho(1:20:end,1:20:end),lat_rho(1:20:end,1:20:end),Windx,Windy,0.5,'k');
    sdate=(datestr(ocean_time(tidx),0));
    date_text=text(-87,41,sdate,'FontSize',20);
    title_t='Significant Wave Height';
    title_text=text(-87,43,title_t,'FontSize',20);
    %hc=colorbar('FontSize',20,'Location','South');
    %set(hc,'Position',[0.8650 0.34 0.0329 0.58])
    caxis([minh maxh+1]);
    dasp(33.5)

    hold on;
    fixpaper2
    drawnow
    %   pause(0.1)
    rootname=(['USEAST_COAWST_',num2str(tidx)]);
    %eval(['imwrite(gcf,''',rootname,'.gif'',''gif'');'])%imwrite will not
    %work when monitor is asleep, writes blank picture
    % Rich Signell's Magick trick (rsignell@usgs.gov) October 6, 2005
    axis off
    pos0=get(0,'screensize');
    posmax=[10 50 pos0(3)-1932 pos0(4)-100];
    % EPSG:4326 is equally spaced lon/lat
    set (gca, 'DataAspectRatio', [1 1 1] );
    set(gcf,'color','w');
    set(gcf,'pos',posmax);
    hold on;
    fixpaper2
    drawnow
    F=getframe(gca);
    pngfile=[rootname '.png'];
    giffile=[rootname '.gif'];
    eval(['print -dpng -r0 ',rootname,'.png;']);

    %          im=F.cdata;
    %          imwrite(im,pngfile,'png');

    ax=axis;
    save ax ax
    %         fid=fopen([rootname '.xml'],'wt');
    %         fprintf(fid,'<BoundingBox>\n');
    %         fprintf(fid,'  <North>\n');
    %         fprintf(fid,'     <Value>%f</Value>\n',ax(4));
    %         fprintf(fid,'  </North>\n');
    %         fprintf(fid,'  <South>\n');
    %         fprintf(fid,'     <Value>%f</Value>\n',ax(3));
    %         fprintf(fid,'  </South>\n');
    %         fprintf(fid,'  <West>\n');
    %         fprintf(fid,'     <Value>%f</Value>\n',ax(1));
    %         fprintf(fid,'  </West>\n');
    %         fprintf(fid,'  <East>\n');
    %         fprintf(fid,'     <Value>%f</Value>\n',ax(2));
    %         fprintf(fid,'  </East>\n');
    %         fprintf(fid,'</BoundingBox>\n');
    %         fclose(fid)

    % write snippet for Google Earth
    %         fid=fopen([rootname '.kml'],'wt');
    %         fprintf(fid,'<LatLonBox>\n');
    %         fprintf(fid,'  <north>%f</north>\n',ax(4));
    %         fprintf(fid,'  <south>%f</south>\n',ax(3));
    %         fprintf(fid,'  <west>%f</west>\n',ax(1));
    %         fprintf(fid,'  <east>%f</east>\n',ax(2));
    %         fprintf(fid,'</LatLonBox>\n');
    %         fclose(fid)
    % convert white to transparent using ImageMagick's "convert" function
    %system(['convert -transparent white ' pngfile ' foo.png']);
    %system(['export PATH=c/programs/ImageMagick-6.4.8/VisualMagick/bin:$PATH$']);
    crop_2_COAWST(pngfile);
    system(['iconvert.exe -transparent white ' pngfile ' foo.gif']);
    movefile('foo.gif',giffile,'f');
    system(['iconvert.exe -transparent white ' pngfile ' foo.png']);
    movefile('foo.png',pngfile,'f');
    frame=frame+1;
end
close all;
ncclose
% end
%crop('D:\Temporary\COAWST\SWAN\run\pngfiles\');