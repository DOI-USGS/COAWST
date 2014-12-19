function [spherical,x,y,bath,Rscope]=read_scope(Gname);

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

% svn $Id: read_scope.m 436 2010-01-02 17:07:55Z arango $
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

[Dnames,Dsizes]=nc_dim(Gname);
ndims=length(Dsizes);
for n=1:ndims,
  dimid=n;
  name=deblank(Dnames(n,:));
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

got.spher =0;  Vname.spher ='spherical';
got.h     =0;  Vname.h     ='h';
got.hraw  =0;  Vname.hraw  ='hraw';
got.rmask =0;  Vname.rmask ='mask_rho';
got.Rscope=0;  Vname.Rscope='scope_rho';
got.rlon  =0;  Vname.rlon  ='lon_rho';
got.rlat  =0;  Vname.rlat  ='lat_rho';
got.xr    =0;  Vname.xr    ='x_rho';
got.yr    =0;  Vname.yr    ='y_rho';

[varnam,nvars]=nc_vname(Gname);
for n=1:nvars,
  name=deblank(varnam(n,:));
  switch name
    case {Vname.spher}
      got.spher=1;
    case {Vname.h}
      got.h=1;
    case {Vname.hraw}
      got.hraw=1;
    case {Vname.rmask}
      got.rmask=1;
    case {Vname.Rscope}
      got.Rscope=1;
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
  spher=nc_read(Gname,Vname.spher);
  if (spher == 'T' | spher == 't'),
    spherical=1;
  end,
end,

% Grid positions at RHO-points.

if (spherical),
  if (got.rlon & got.rlat),
    x=nc_read(Gname,Vname.rlon);
    y=nc_read(Gname,Vname.rlat);
  else,
    [y,x]=meshgrid(1:Jm,1:Im);
  end,
else,
  if (got.xr & got.yr),
    x=nc_read(Gname,Vname.xr);
    y=nc_read(Gname,Vname.yr);
  else,
    [y,x]=meshgrid(1:Jm,1:Im);
  end,
end,

% Scope mask on RHO-points.

if (got.Rscope),
  Rscope=nc_read(Gname,Vname.Rscope);
else,
  if (got.rmask),
    Rscope=nc_read(Gname,Vname.rmask);
    size(Rscope)
  else,
    Rscope=ones(size(x));
  end,
end,

% Bathymetry.

if (got.hraw),
  bath=nc_read(Gname,Vname.hraw,1);
else,
  if (got.h),
    bath=nc_read(Gname,Vname.h);
  else,
    bath=zeros(size(x));
  end,
end,

return
