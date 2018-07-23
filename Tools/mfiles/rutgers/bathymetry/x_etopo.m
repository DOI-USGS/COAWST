function [lon,lat,h]=x_etopo(Llon, Rlon, Blat, Tlat, Fname);

%
% X_ETOPO:  Extract requested bathymetry from ETOPO dataset
%
% [lon,lat,h]=x_etopo5(Llon, Rlon, Blat, Tlat, Fname)
%
% This function extract bathymetry data from a ETOPO dataset NetCDF file
% in the specified geographical area.
%
%
% On Input:
%
%    Llon          Left   corner longitude (degrees_east)
%    Rlon          Right  corner longitude (degrees_east)
%    Blat          Bottom corner latitude (degrees_north)
%    Tlat          Top    corner latitude (degrees_north)
%
% On Output:
%
%    h             Bathymetry (m, matrix)
%    lon           Longitude (degree_east, matrix)
%    lat           Latitude (degree_north, matrix)
%

% svn $Id: x_etopo.m 895 2018-02-11 23:15:37Z arango $
%=========================================================================%
%  Copyright (c) 2002-2018 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  Read in ETOPO data longitude and latitude.  The ETOPO dataset longitude
%  in the provided file goes from -180 to 180.  
%--------------------------------------------------------------------------

x=nc_read(Fname,'topo_lon');
y=nc_read(Fname,'topo_lat');

%--------------------------------------------------------------------------
%  If geographical area is a both sides of the day line (Lon=-180=180),
%  we need to convert to east longitudes: 0 to 360.
%--------------------------------------------------------------------------

if (abs(Llon) > 180 || abs(Rlon) > 180),
  ind=find(x < 0);
  x(ind)=360+x(ind);
end

%--------------------------------------------------------------------------
%  Determine indices to extract.
%--------------------------------------------------------------------------

I=find(x >= Llon & x <= Rlon);
X=x(I);
[X,Isort]=sort(X);
I=I(Isort);

J=find(y >= Blat & y <= Tlat);
Jmin=min(J);
Jmax=max(J);
Y=y(J); Y=Y';

%--------------------------------------------------------------------------
%  Read in bathymetry.
%--------------------------------------------------------------------------

topo=nc_read(Fname,'topo');
h=topo(I,Jmin:Jmax);
[Im,Jm]=size(h);

%--------------------------------------------------------------------------
%  Set coordinates of extracted bathymetry.
%--------------------------------------------------------------------------

lon=repmat(X,[1 Jm]);
lat=repmat(Y,[Im 1]);

return
