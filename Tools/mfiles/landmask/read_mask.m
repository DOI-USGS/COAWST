function [spherical,x,y,bath,rmask]=read_mask(Gname);

%
% READ_SCOPE:  Read ROMS Land/Sea mask
%
% [spherical,x,y,bath,rmask]=read_mask(Gname)
%
% This routine reads in domain grid, bathymetry, and Land/Sea mask at
% from GRID NetCDF file.
%
% On Input:
%
%    Gname       GRID NetCDF file name (character string).
%
% On Output:
%
%    spherical   grid type switch (logical):
%                  spherical=1, spherical grid set-up.
%                  spherical=0, Cartesian grid set-up.
%    x           X-location of RHO-points (real matrix).
%    y           Y-location of RHO-points (real matrix).
%    bath        raw bathymetry at RHO-points (real matrix; meters).
%    rmask       Land/Sea mask on RHO-points (real matrix):
%                  rmask=0 land, rmask=1 Sea.
%

% svn $Id: read_mask.m 436 2010-01-02 17:07:55Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%---------------------------------------------------------------------------
% Inquire about spatial dimensions.
%---------------------------------------------------------------------------

Dname.xr='xi_rho';
Dname.yr='eta_rho';

%[Dnames,Dsizes]=nc_dim(Gname);
finfo=ncinfo(Gname);
Dnames={finfo.Dimensions.Name};
Dsizes=[finfo.Dimensions.Length];
ndims=length(Dsizes);
for n=1:ndims,
  dimid=n;
  name=deblank(char(Dnames(n)));
  switch name
    case {Dname.xr}
      Im=Dsizes(n);
    case {Dname.yr}
      Jm=Dsizes(n);
  end,
end,

%---------------------------------------------------------------------------
% Inquire grid NetCDF file about mask variables.
%---------------------------------------------------------------------------

got.spher=0;  Vname.spher='spherical';
got.h    =0;  Vname.h    ='h';
got.hraw =0;  Vname.hraw ='hraw';
got.rmask=0;  Vname.rmask='mask_rho';
got.rlon =0;  Vname.rlon ='lon_rho';
got.rlat =0;  Vname.rlat ='lat_rho';
got.xr   =0;  Vname.xr   ='x_rho';
got.yr   =0;  Vname.yr   ='y_rho';

%[varnam,nvars]=nc_vname(Gname);
Dnames={finfo.Variables.Name};
nvars=length(Dnames);

for n=1:nvars,
  name=deblank(char(Dnames(n)));
  switch name
    case {Vname.spher}
      got.spher=1;
    case {Vname.h}
      got.h=1;
%    case {Vname.hraw}
%      got.hraw=1;
    case {Vname.rmask}
      got.rmask=1;
    case {Vname.xr}
      got.xr=1;
    case {Vname.yr}
      got.yr=1;
  end,
end,

%---------------------------------------------------------------------------
% Read in relevant Land/Sea mask variables.
%---------------------------------------------------------------------------

% Spherical switch.

spherical=0;
if (got.spher),
  spher=ncread(Gname,Vname.spher);
  if (spher == 'T' | spher == 't'),
    spherical=1;
  end,
end,

% Grid positions at RHO-points.

if (spherical),
  if (got.rlon & got.rlat),
    x=ncread(Gname,Vname.rlon);
    y=ncread(Gname,Vname.rlat);
  else,
    [y,x]=meshgrid(1:Jm,1:Im);
  end,
else,
  if (got.xr & got.yr),
    x=ncread(Gname,Vname.xr);
    y=ncread(Gname,Vname.yr);
  else,
    [y,x]=meshgrid(1:Jm,1:Im);
  end,
end,

% Mask on RHO-points.

if (got.rmask),
  rmask=ncread(Gname,Vname.rmask);
else,
  rmask=ones(size(x));
end,

% Bathymetry.

if (got.hraw),
  bath=ncread(Gname,Vname.hraw,1);
else,
  if (got.h),
    bath=ncread(Gname,Vname.h);
  else,
    bath=zeros(size(x));
  end,
end,

return
