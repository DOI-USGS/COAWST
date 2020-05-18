function [spherical,x,y,bath,rmask]=read_mask(ncfile)

%
% READ_MASK:  Read ROMS Land/Sea mask
%
% [spherical,x,y,bath,rmask]=read_mask(ncfile)
%
% This routine reads in domain grid, bathymetry, and Land/Sea mask at
% from GRID NetCDF file.
%
% On Input:
%
%    ncfile      GRID NetCDF file name (string)
%
% On Output:
%
%    spherical   Spherical switch (logical):
%                  spherical = true,  spherical grid set-up
%                  spherical = false, Cartesian grid set-up
%    x           X-location of RHO-points (2D array)
%    y           Y-location of RHO-points (2D array)
%    bath        raw bathymetry at RHO-points (2D array; meters)
%    rmask       Land/Sea mask on RHO-points (2D array):
%                  rmask=0 land, rmask=1 Sea.
%

% svn $Id: read_mask.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%-------------------------------------------------------------------------
% Inquire about spatial dimensions.
%-------------------------------------------------------------------------

Dname.xr='xi_rho';
Dname.yr='eta_rho';

D=nc_dinfo(ncfile);
ndims=length(D);
for n=1:ndims,
  name=char(D(n).Name);
  switch name
    case {Dname.xr}
      Im=D(n).Length;
    case {Dname.yr}
      Jm=D(n).Length;
  end
end

%-------------------------------------------------------------------------
% Inquire grid NetCDF file about mask variables.
%-------------------------------------------------------------------------

got.spher=false;  Vname.spher='spherical';
got.h    =false;  Vname.h    ='h';
got.hraw =false;  Vname.hraw ='hraw';
got.rmask=false;  Vname.rmask='mask_rho';
got.rlon =false;  Vname.rlon ='lon_rho';
got.rlat =false;  Vname.rlat ='lat_rho';
got.xr   =false;  Vname.xr   ='x_rho';
got.yr   =false;  Vname.yr   ='y_rho';

V=nc_vnames(ncfile);
nvars=length(V.Variables);
for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch name
    case {Vname.spher}
      got.spher=true;
    case {Vname.h}
      got.h=true;
    case {Vname.hraw}
      if ~any(V.Variables(n).Size == 0),
        got.hraw=true;
      end
    case {Vname.rmask}
      got.rmask=true;
    case {Vname.xr}
      got.xr=true;
    case {Vname.yr}
      got.yr=true;
  end
end

%-------------------------------------------------------------------------
% Read in relevant Land/Sea mask variables.
%-------------------------------------------------------------------------

% Spherical switch.

if (got.spher),
  spherical=nc_read(ncfile,Vname.spher);
  if (ischar(spherical)),
    if (spherical == 'T' || spherical == 't')
      spherical = true;
    else
      spherical = false;
    end
  end
end

% Grid positions at RHO-points.

if (spherical),
  if (got.rlon && got.rlat),
    x=nc_read(ncfile,Vname.rlon);
    y=nc_read(ncfile,Vname.rlat);
  else
    [y,x]=meshgrid(1:Jm,1:Im);
  end
else
  if (got.xr && got.yr),
    x=nc_read(ncfile,Vname.xr);
    y=nc_read(ncfile,Vname.yr);
  else
    [y,x]=meshgrid(1:Jm,1:Im);
  end
end

% Mask on RHO-points.

if (got.rmask),
  rmask=nc_read(ncfile,Vname.rmask);
else
  rmask=ones(size(x));
end

% Bathymetry.

if (got.hraw),
  bath=nc_read(ncfile,Vname.hraw,1);
else
  if (got.h),
    bath=nc_read(ncfile,Vname.h);
  else
    bath=zeros(size(x));
  end
end

return
