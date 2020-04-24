%
%  EXTRACT_BATH:  Driver script to extract bathymetry data
%
%  This user modifiable script can extract bathymetry from ETOPO-5 database
%  at the specified coordinates.
%

% svn $Id: extract_bath.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%job='seagrid';          %  prepare bathymetry for seagrid
 job='netcdf';           %  Extract a bathymetry NetCDF file

database='etopo';

switch job,
  case 'seagrid'
    Oname='uswest_bath.mat';
  
  case 'netcdf'
    Oname='uswest_bath.nc';
end,

switch database,
  case 'etopo'
    Bname='etopo5.nc';
end,

Llon=-134.0;              % Left   corner longitude     % US west Coast
Rlon=-118.0;              % Right  corner longitude
Blat=33.0;                % Bottom corner latitude
Tlat=49.0;                % Top    corner latitude

%-----------------------------------------------------------------------
%  Read and extract bathymetry.
%-----------------------------------------------------------------------

switch database,
  case 'etopo'
    [lon,lat,h]=x_etopo(Llon,Rlon,Blat,Tlat,Bname);
end,

%-----------------------------------------------------------------------
%  Process extracted data for requested task.
%-----------------------------------------------------------------------

switch job,
  case 'seagrid'
    xbathy=reshape(lon,1,prod(size(lon)))';
    ybathy=reshape(lat,1,prod(size(lat)))';
    zbathy=reshape(h,1,prod(size(h)))';
    zbathy=-zbathy;
    ind=find(zbathy<0);
    if (~isempty(ind)),
      zbathy(ind)=0;
    end;
    save(Oname,'xbathy','ybathy','zbathy');
  
  case 'netcdf'
    [Im,Jm]=size(h);  
    status=c_bath(Im,Jm,Oname);
    status=nc_write(Oname,'spherical',1);
    status=nc_write(Oname,'lon',lon);
    status=nc_write(Oname,'lat',lat);
    status=nc_write(Oname,'hraw',h,1);
end,
