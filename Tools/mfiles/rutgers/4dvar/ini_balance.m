function [K]=ini_balance(Gname,Hname,TimeRec);

%
% INI_BALANCE:  Read and initialize variables used in the balance operator
%
% [K]=ini_balance(Gname,Hname)
%
% This function computes the balanced, baroclinic U- and V-momentum 
% anomalies using geostrophic balance for the error covariance balance
% operator, K^(-1). It uses ROMS pressure gradient "prsgrd31.h" formulation
% to compute the dynamic pressure from the balance density anomaly computed
% in function "rho_balance".
%
% On Input:
%
%    Gname     ROMS Grid NetCdF file name (character string)
%    Hname     ROMS History NetCDF file name (character string)
%    TimeRec   Time record of History NetCDF file to use (integer)
%
% On Output:
%
%    K         Balance operator structure array:
%
%              K.Lm           Number of interior points in XI-direction
%              K.Mm           Number of interior points in ETA-direction
%              K.N            Number of vertical level
%              K.Niter        Number of iterations in the biconjugate
%                               gradient solver.
%              K.EW_periodic  Switch for East-West periodic boundaries
%              K.NS_periodic  Switch for North-South periodic boundaries
%              K.g            Accelerarion due to gravity (m/s2)
%              K.rho0         Mean density (Kg/m3) used when the Boussinesq
%                               approximation is inferred
%              K.ml_depth     Mixed-layer depth (100 m)
%              K.dTdz_min     Minimum dT/dz allowed (0.001 Celsius/m)
%              K.Norder       Shapiro filter order to use
%              K.Fcoef        Shapiro filter coefficients
%              K.h            Bathymetry (m) at RHO-points
%              K.f            Coriolis parameter (1/s) at RHO-points
%              K.pm           Curvilinear X-coordinate metric (1/m)
%              K.pn           Curvilinear Y-coordinate metric (1/m)
%              K.pmon_u       Curvilinear metric pm/pn at U-points
%              K.pnom_v       Curvilinear metric pn/pm at V-points
%              K.Pmask        Land/Sea mask on PSI-points
%              K.rmask        Land/Sea mask on RHO-points
%              K.umask        Land/Sea mask on U-points
%              K.vmask        Land/Sea mask on V-points
%              K.alpha        Surface thermal expansion coefficient
%                               (1/Celsius)
%              K.beta         Surface saline contraction coefficient
%                               (nondimensional)
%              K.rho          Basic state in situ density (kg/m3)
%              K.Hz           Vertical level thicknesses (m)
%              K.Zr           Depths at vertical RHO-points (m, negative)
%              K.Zw           Depths at vertical W-points   (m, negative)
%

% svn $Id: ini_balance.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Load input NetCDF file names.

K.Gname=Gname;
K.Hname=Hname;

%  Set number of grid points.

[Lp,Mp]=size(nc_read(Gname,'h'));
L=Lp-1; Lm=L-1;
M=Mp-1; Mm=M-1;

K.Lm=Lm;
K.Mm=Mm;

[N]=length(nc_read(Hname,'s_rho'));

K.N=N;

%  Load History record processed.

K.TimeRec=TimeRec;

%  Set number of iterations in the biconjugate gradient algorithm for the
%  SSH elliptic equation solution.

K.Niter=200;

%  Set switches for periodic boundary conditions. Inquire History NetCDF
%  global attribute "CPP_options" to see if either EW_PERIODIC or
%  NS_PERIODIC is defined.

[CPP]=nc_getatt(Hname,'CPP_options');

K.EW_periodic=0;
K.NS_periodic=0;

if (~isempty(CPP)),
 foundit=strfind(CPP,'EW_PERIODIC');
 if isempty(~foundit),
   K.EW_PERIODIC=1;
 end,

 foundit=strmatch(CPP,'NS_PERIODIC');
 if isempty(~foundit),
   K.NS_PERIODIC=1;
 end,
end,

%  Set accelerarion due to gravity (m/s2).

K.g=9.81;

%  Read mean density (Kg/m3) used when the Boussinesq approximation is
%  inferred.

K.rho0=nc_read(Hname,'rho0');

%  Set default mixed-layer depth (m, positive) used to compute balance
%  salinity from temperature.

K.ml_depth=100;

%  Set default minimum dT/dz value (Celsius/m) to consider when computing
%  balance salinity from temperature.

K.dTdz_min=0.001;

%  Set Shapiro filter coeficients and order to use.

K.Norder=2;

K.Fcoef=[2.500000D-1   6.250000D-2    1.562500D-2    3.906250D-3     ...
         9.765625D-4   2.44140625D-4  6.103515625D-5 1.5258789063D-5 ...
         3.814697D-6   9.536743D-7    2.384186D-7    5.960464D-8     ...
         1.490116D-8   3.725290D-9    9.313226D-10   2.328306D-10    ...
         5.820766D-11  1.455192D-11   3.637979D-12   9.094947D-13];

%  Read in bathymetry (m) at RHO-points.

K.h=nc_read(Gname,'h');

%  Read in Coriolis parameter (1/s) at RHO-points.

K.f=nc_read(Gname,'f');

%  Read in curvilinear metrics (1/m).

K.pm=nc_read(Gname,'pm');
K.pn=nc_read(Gname,'pn');

%  Compute pm/pn and pn/pm metrics (nondimensional)

K.pmon_u=(K.pm(1:L,1:Mp)+K.pm(2:Lp,1:Mp))./                          ...
         (K.pn(1:L,1:Mp)+K.pn(2:Lp,1:Mp));

K.pnom_v=(K.pn(1:Lp,1:M)+K.pn(1:Lp,2:Mp))./                          ...
         (K.pm(1:Lp,1:M)+K.pm(1:Lp,2:Mp));
	      
%  Read in or set land/sea masking.

V=nc_vnames(Gname);
nvars=length(V.Variables);

for n=1:nvars
  name=char(V.Variables(n).Name);
  switch name
    case 'mask_psi'
      K.pmask=nc_read(Gname,'mask_psi');
    case 'mask_rho'
      K.rmask=nc_read(Gname,'mask_rho');
    case 'mask_u'
      K.umask=nc_read(Gname,'mask_u');
    case 'mask_v'
      K.vmask=nc_read(Gname,'mask_v');
  end,
end,

if (~isfield(K,'pmask')),
  K.pmask=ones(size([L,M]));
end,

if (~isfield(K,'rmask')),
  K.rmask=ones(size([Lp,Mp]));
end,

if (~isfield(K,'umask')),
  K.umask=ones(size([L,Mp]));
end,

if (~isfield(K,'vmask')),
  K.vmask=ones(size([Lp,M]));
end,

%  Compute surface thermal expansion (1/Celsius), saline contraction 
%  coefficients, and in situ density (kg/m3).

F=eos(Gname,Hname,TimeRec);

K.alpha=squeeze(F.alpha(:,:,N)).*K.rmask;
K.beta=squeeze(F.beta(:,:,N)).*K.rmask;
K.rho=(F.den-1000).*repmat(K.rmask,[1 1 N]);

%  Compute vertical depths (m, negative).  Notice that "tindex" is set
%  to zero for zero free-surface since we want the rest state depths.

tindex=0;

K.Zr=depths(Hname, Gname, 1, 0, tindex);
K.Zw=depths(Hname, Gname, 5, 0, tindex);

%  Compute vertical level thicknesses (m, positive).

K.Hz(1:Lp,1:Mp,1:N)=K.Zw(1:Lp,1:Mp,2:N+1)-K.Zw(1:Lp,1:Mp,1:N);

return
