function rivers_latlon_to_IJ(Gname,Rname)
river_coord=load(Rname);
rivers_lon=river_coord(:,1);
rivers_lat=river_coord(:,2); 
rlon=nc_read(Gname,'lon_rho');
rlat=nc_read(Gname,'lat_rho');
[Lp,Mp]=size(rlon);
[y,x]=meshgrid(1:Mp,1:Lp);
LonAvg=mean(rlon(:));
LatAvg=mean(rlat(:));
I= griddata(rlon-LonAvg,rlat-LatAvg,x,rivers_lon-LonAvg,...
    rivers_lat-LatAvg,'linear');
J= griddata(rlon-LonAvg,rlat-LatAvg,y,rivers_lon-LonAvg,...
    rivers_lat-LatAvg,'linear');
I=int16(I-1);
J=int16(J-1);
IJ= double([I J]);
save('rivers_IJ.txt', 'IJ','-ASCII');
type rivers_IJ.txt;
end