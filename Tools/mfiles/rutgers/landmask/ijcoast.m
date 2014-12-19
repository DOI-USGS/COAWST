function [C]=ijcoast(Gname, Cname);

%
% IJCOAST:  Converts coastline (lon,lat) to fractional coordinates
%
% [C]=ijcoast(Gname,Cname)
%
% This script converts coastline data to GRID indices for a particular
% application.  This will used in the Land/Sea masking tools.
%
% On Input:
%
%    Gname       NetCDF file name (character string).
%    Cname       Coastline file name (character string).
%
% On Output:
%
%    C           Coastline indices (structure array):
%                  C.grid    => Grid file name.
%                  C.coast   => Coastline file name.
%                  C.indices => Coastline indices file name.
%                  C.lon     => Coastline longitudes.
%                  C.lat     => Coastline latitudes.
%                  C.Icst    => Coastline I-grid coordinates, (0:L).
%                  C.Jcst    => Coastline J-grid coordinates, (0:M).
%

% svn $Id: ijcoast.m 436 2010-01-02 17:07:55Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

EXTRACT=0;
method='linear';

%---------------------------------------------------------------------------
% Inquire NetCDF about coastline data.
%---------------------------------------------------------------------------

got_coast=0;
if (nargin < 2),
% [varnam,nvars]=nc_vname(Gname);
  finfo=ncinfo(grid_file);
  Dnames={finfo.Variables.Name};
  got_coast=strmatch('lon_coast',Dnames,'exact');
  got_coast=got_coast+strmatch('lat_coast',Dnames,'exact');
end,

C.grid=Gname;
if (~got_coast),
  C.coast=Cname;
end,
C.indices='ijcoast.mat';

%---------------------------------------------------------------------------
% Read in grid coordinates at rho-points.
%---------------------------------------------------------------------------

rlon=ncread(Gname,'lon_rho');
rlat=ncread(Gname,'lat_rho');

[Lp,Mp]=size(rlon);
L=Lp-1;
M=Mp-1;

%---------------------------------------------------------------------------
% Read in coastline data.
%---------------------------------------------------------------------------

if (got_coast),
  Clon=ncread(Gname,'lon_coast');
  Clat=ncread(Gname,'lat_coast');
else
  load(Cname);
  Clon=lon;
  Clat=lat;
end,

%---------------------------------------------------------------------------
% Extract need coasline data.
%---------------------------------------------------------------------------

if (EXTRACT),

  dx=5*abs(mean(mean(gradient(rlon))));
  dy=5*abs(mean(mean(gradient(rlat))));

  C.Llon=min(min(rlon));  C.Llon=C.Llon-dx;
  C.Rlon=max(max(rlon));  C.Rlon=C.Rlon+dx;
  C.Blat=min(min(rlat));  C.Blat=C.Blat-dy;
  C.Tlat=max(max(rlat));  C.Tlat=C.Tlat+dy;

  ind=find(Clon > C.Rlon | Clon < C.Llon | Clat > C.Tlat | Clat < C.Blat);
  clon=Clon;
  clat=Clat;
  if (~isempty(ind)),
    clon(ind)=[];
    clat(ind)=[];
  end,

  C.lon=[NaN; clon; NaN];
  C.lat=[NaN; clat; NaN];

else,

  clon=Clon;
  clat=Clat;

  C.lon=clon;
  C.lat=clat;

end,

clear Clon Clat clon clat

%---------------------------------------------------------------------------
% Interpolate coastline to grid units. Notice that the mean grid
% longitude and latitude values are substracted to avoid roundoff
% errors from triangulation.
%---------------------------------------------------------------------------

disp(['Converting coastline (lon,lat) to (I,J) grid indices.']);

[y,x]=meshgrid(1:Mp,1:Lp);

LonAvg=mean(rlon(:));
LatAvg=mean(rlat(:));

C.Icst=griddata(rlon-LonAvg,rlat-LatAvg,x,C.lon-LonAvg,C.lat-LatAvg,method);
C.Jcst=griddata(rlon-LonAvg,rlat-LatAvg,y,C.lon-LonAvg,C.lat-LatAvg,method);

%  Substrat one to have indices in the range (0:L,0:M).

C.Icst=C.Icst-1;
C.Jcst=C.Jcst-1;

%---------------------------------------------------------------------------
%  Save sctructure array into a Matlab file.
%---------------------------------------------------------------------------

disp(['Saving coastline (I,J) into file: ',C.indices]);

save(C.indices,'C');

return
