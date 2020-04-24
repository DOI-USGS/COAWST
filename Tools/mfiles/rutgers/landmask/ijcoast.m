function [C]=ijcoast(ncfile, coast_file)

%
% IJCOAST:  Converts coastline (lon,lat) to fractional coordinates
%
% [C]=ijcoast(ncfile,coast_file)
%
% Converts coastline data to ROMS (I,J) fractional coordinates to
% facilitate Land/Sea masking editing with 'editmask'.
%
% On Input:
%
%    ncfile      NetCDF file name (string)
%    coast_file  Coastline file name (string)
%
% On Output:
%
%    C           Coastline indices (structure array):
%                  C.grid    => Grid file name
%                  C.coast   => Coastline file name
%                  C.indices => Coastline indices file name
%                  C.lon     => Coastline longitudes
%                  C.lat     => Coastline latitudes
%                  C.Icst    => Coastline I-grid coordinates, (0:L)
%                  C.Jcst    => Coastline J-grid coordinates, (0:M)
%

% svn $Id: ijcoast.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

EXTRACT=false;
method='linear';

%-------------------------------------------------------------------------
% Inquire NetCDF about coastline data.
%-------------------------------------------------------------------------

got_coast=false;
if (nargin < 2)
  S=nc_vnames(ncfile);
  vnames={S.Variables.Name};
  got_Clon=any(strcmp(vnames,'lon_coast'));
  got_Clat=any(strcmp(vnames,'lat_coast'));
  got_coast=got_Clon && got_Clat;
end

C.grid=ncfile;
if (~got_coast)
  C.coast=coast_file;
end
C.indices='ijcoast.mat';

%--------------------------------------------------------------------------
% Read in grid coordinates at rho-points.
%--------------------------------------------------------------------------

rlon=nc_read(ncfile,'lon_rho');
rlat=nc_read(ncfile,'lat_rho');

[Lp,Mp]=size(rlon);

%--------------------------------------------------------------------------
% Read in coastline data.
%--------------------------------------------------------------------------

if (got_coast)
  Clon=nc_read(ncfile,'lon_coast');
  Clat=nc_read(ncfile,'lat_coast');
else
  load(coast_file);
  Clon=lon;
  Clat=lat;
end

%--------------------------------------------------------------------------
% Extract need coasline data.
%--------------------------------------------------------------------------

if (EXTRACT)

  dx=5*abs(mean(mean(gradient(rlon))));
  dy=5*abs(mean(mean(gradient(rlat))));

  C.Llon=min(min(rlon));  C.Llon=C.Llon-dx;
  C.Rlon=max(max(rlon));  C.Rlon=C.Rlon+dx;
  C.Blat=min(min(rlat));  C.Blat=C.Blat-dy;
  C.Tlat=max(max(rlat));  C.Tlat=C.Tlat+dy;

  ind=find(Clon > C.Rlon | Clon < C.Llon | Clat > C.Tlat | Clat < C.Blat);
  clon=Clon;
  clat=Clat;
  if (~isempty(ind))
    clon(ind)=[];
    clat(ind)=[];
  end

  C.lon=[NaN; clon; NaN];
  C.lat=[NaN; clat; NaN];

else

  clon=Clon;
  clat=Clat;

  C.lon=clon;
  C.lat=clat;

end

clear Clon Clat clon clat

%--------------------------------------------------------------------------
% Interpolate coastline to grid units. Notice that the mean grid
% longitude and latitude values are substracted to avoid roundoff
% errors from triangulation.
%--------------------------------------------------------------------------

disp('Converting coastline (lon,lat) to (I,J) grid indices.');

[y,x]=meshgrid(1:Mp,1:Lp);

C.Icst=griddata(rlon,rlat,x,C.lon,C.lat,method);
C.Jcst=griddata(rlon,rlat,y,C.lon,C.lat,method);

%  Substrat one to have indices in the range (0:L,0:M).

C.Icst=C.Icst-1;
C.Jcst=C.Jcst-1;

%--------------------------------------------------------------------------
%  Save sctructure array into a Matlab file.
%--------------------------------------------------------------------------

disp(['Saving coastline (I,J) into file: ',C.indices]);

save(C.indices,'C');

return
