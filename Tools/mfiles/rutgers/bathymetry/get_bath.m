function [lon,lat,h]=get_bath(Llon,Rlon,Blat,Tlat,varargin)
  
%
% GET_BATH:  Extract bathymetry and elevation from a ETOPO NetCDF file
%
% [lon,lat,h]=get_bath(Llon,Rlon,Blat,Tlat,OutFile,InpFile)
%
% This function extracts bathymetry and elevation data from a ETOPO
% NetCDF filw for the box corners bounded by (Llon,Blat) and (Rlon,Tlat).
%
%                 ______ (Rlon, Tlat)
%                |      |
%                |      |
%                |______|
%    (Llon, Blat)    
%
% On Input:
%
%    Llon         Box left-edge   longitude (degrees, -180 - 180)
%
%    Rlon         Box right-edge  longitude (degrees, -180 - 180)
%
%    Blat         Box bottom-edge latitude  (degress,  -90 - 90 )
%
%    Tlat         Box top-edge    latitude  (degress,  -90 - 90 )
%
%    OutFile      Output extracted bathymetry file (Optional, string). The
%                   following values are possible:
%
%                   OutFile = []          no output file (default)
%                   OutFile = 'xxxx.mat'  SeaGrid Matlab file
%                   OutFile = 'xxxx.nc'   NetCDF file
%
%    InpFile      Input ETOPO file (Optional, string)
%
%                   InpFile = '~/ocean/repository/matlab/bathymetry/etopo5.nc'
%                             (default)
%
% On Ouput:
%
%    lon          Extracted coastline longitude       
%
%    lat          Extracted coastline latitude       
%
%    h            Extracted bathymetry and elevation (m)       
%
 
% svn $Id: get_bath.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Set optional arguments

OutFile = [];
InpFile = '~/ocean/repository/matlab/bathymetry/etopo5.nc';

switch numel(varargin)
 case 1
   OutFile = varargin{1};
 case 2
   OutFile = varargin{1};
   InpFile = varargin{2};
end

%--------------------------------------------------------------------------
% Extract coastlines from GSHHS database.
%--------------------------------------------------------------------------

[lon,lat,h]=x_etopo(Llon,Rlon,Blat,Tlat,InpFile);

%---------------------------------------------------------------------------
%  Save extrated bathymetry.
%---------------------------------------------------------------------------

if (~isempty(OutFile)),
  [~, ~, extension] = fileparts(OutFile);
  switch extension,
    case '.mat'                                          % seagrid
      xbathy = reshape(lon,1,prod(size(lon)))';
      ybathy = reshape(lat,1,prod(size(lat)))';
      zbathy = reshape(h,1,prod(size(h)))';
      zbathy = -zbathy;
      ind = find(zbathy<0);
      if (~isempty(ind)),
        zbathy(ind) = 0;
      end;
      save (OutFile, 'xbathy', 'ybathy', 'zbathy');
    case '.nc'                                           % NetCDF
      [Im,Jm] = size(h);
      status = c_bath(Im, Jm, OutFile);
      status = nc_write(OutFile, 'spherical', 1);
      status = nc_write(OutFile, 'lon', lon);
      status = nc_write(OutFile, 'lat', lat);
      status = nc_write(OutFile, 'hraw', h, 1);
  end
end

return
