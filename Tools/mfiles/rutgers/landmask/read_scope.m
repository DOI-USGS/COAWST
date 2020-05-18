function [spherical,x,y,bath,Rscope]=read_scope(Gname)

%
% READ_SCOPE:  Read ROMS Adjoint sensitivity scope mask
%
% [spherical,x,y,bath,Rscope]=read_scope(Gname)
%
% This routine reads in domain grid, bathymetry, and adjoint
% sensitivity scope mask from GRID NetCDF file.
%
% On Input:
%
%    Gname       GRID NetCDF file name (character string).
%
% On Output:
%
%    spherical   Grid type switch (logical):
%                  spherical=1, spherical grid set-up.
%                  spherical=0, Cartesian grid set-up.
%    x           X-location of RHO-points (real matrix).
%    y           Y-location of RHO-points (real matrix).
%    bath        Raw bathymetry at RHO-points (real matrix; meters).
%    Rscope      Adjoint sensitivity scope mask on RHO-points (real
%                  matrix): Rscope=0 inactive, Rscope=1 active.
%

% svn $Id: read_scope.m 996 2020-01-10 04:28:56Z arango $
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

D=nc_dinfo(Gname);
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

got.spher =false;  Vname.spher ='spherical';
got.h     =false;  Vname.h     ='h';
got.hraw  =false;  Vname.hraw  ='hraw';
got.rmask =false;  Vname.rmask ='mask_rho';
got.Rscope=false;  Vname.Rscope='scope_rho';
got.rlon  =false;  Vname.rlon  ='lon_rho';
got.rlat  =false;  Vname.rlat  ='lat_rho';
got.xr    =false;  Vname.xr    ='x_rho';
got.yr    =false;  Vname.yr    ='y_rho';

V=nc_vnames(Gname);
nvars=length(V.Variables);
for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch name
    case {Vname.spher}
      got.spher=true;
    case {Vname.h}
      got.h=true;
    case {Vname.hraw}
      got.hraw=true;
    case {Vname.rmask}
      got.rmask=true;
    case {Vname.Rscope}
      got.Rscope=true;
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

spherical=false;
if (got.spher),
  spher=nc_read(Gname,Vname.spher);
  if (spher == 'T' || spher == 't'),
    spherical=true;
  end
end

% Grid positions at RHO-points.

if (spherical),
  if (got.rlon && got.rlat),
    x=nc_read(Gname,Vname.rlon);
    y=nc_read(Gname,Vname.rlat);
  else
    [y,x]=meshgrid(1:Jm,1:Im);
  end
else
  if (got.xr && got.yr),
    x=nc_read(Gname,Vname.xr);
    y=nc_read(Gname,Vname.yr);
  else
    [y,x]=meshgrid(1:Jm,1:Im);
  end
end

% Scope mask on RHO-points.

if (got.Rscope),
  Rscope=nc_read(Gname,Vname.Rscope);
else
  if (got.rmask),
    Rscope=nc_read(Gname,Vname.rmask);
  else
    Rscope=ones(size(x));
  end
end

% Bathymetry.

if (got.hraw),
  bath=nc_read(Gname,Vname.hraw,1);
else
  if (got.h),
    bath=nc_read(Gname,Vname.h);
  else
    bath=zeros(size(x));
  end
end

return
