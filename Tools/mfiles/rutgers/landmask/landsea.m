function [F]=landsea(ncfile, varargin)

% LANDSEA:  Computes and writes ROMS Land/Sea masks
%
% [F]=landsea(ncfile)
% [F]=landsea(ncfile, database)
% [F]=landsea(ncfile, clon, clat)
%
% Computes automatically the Land/Sea mask for an application grid and
% writes it into input NetCDF file. The Land/Sea mask is computed from
% closed coastlines polygoms extracted from the GSHHS dataset.
%
% On Input:
%
%    ncfile      NetCDF file name (string)
%
%    database    GSHHS database (character, OPTIONAL)
%                  'f'    full resolution
%                  'h'    high resolution
%                  'i'    intermediate resolution (default)
%                  'l'    load resolution
%                  'c'    crude resolution
% or
%    clon        Coastline longitude (real vector, OPTIONAL)
%    clat        Coastline latitude  (real vector, OPTIONAL)
%
% On Output:
%
%    F           Land/Sea mask (structure array):
%                  F.mask_rho =>  Mask at rho-points
%                  F.mask_psi =>  Mask at psi-points
%                  F.mask_u   =>  Mask at u-points
%                  F.mask_v   =>  Mask at v-points
%
% Calls:   MEXNC (Interface to NetCDF library using Matlab):
%          nc_dinfo, nc_read, nc_vdef, nc_vnames, nc_write
%
%          r_gshhs, x_gshhs
%
% The masking algorithm is adapted from Charles R. Denham routine
% "domask.m" from the SeaGrid package.
%

% svn $Id: landsea.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

Tstart=tic;                         % time profile

extract_coast=false;                % extract coastlines from GSHHS switch

Cfile=which('gshhs_h.b','-ALL');    % select first directory found
DIR=fileparts(Cfile{1});            % others are shadowed

F.ncfile=ncfile;

switch numel(varargin)
  case 0
    extract_coast=true;
    name='gshhs_i.b';               % intermediate resolution
    Cname=fullfile(DIR,name);
    F.Cname=Cname;
  case 1
    extract_coast=true;
    database=varargin{1};
    switch database
      case 'f'                      % full resolution
        name='gshhs_f.b';
      case 'h'                      % high resolution
        name='gshhs_h.b';
      case 'i'                      % intermediate resolution
        name='gshhs_i.b';
      case 'l'                      % low resolution
        name='gshhs_l.b';
      case 'c'                      % crude resolution
        name='gshhs_c.b';
      otherwise
        error(['LANDSEA: illegal GSHHS dataset resolution, ',database])
    end
    Cname=fullfile(DIR,name);
    F.Cname=Cname;
 case 2
   clon=varargin{1};
   clat=varargin{2};
end

% Open GRID NetCDF file for writing.

ncid=netcdf.open(ncfile, 'WRITE');

%-----------------------------------------------------------------------
% Inquire grid NetCDF file about mask variables.
%-----------------------------------------------------------------------

VarList={'mask_rho', 'mask_psi', 'mask_u', 'mask_v'};

for var = VarList
  field=char(var);
  got.(field)=false;
  Vname.(field)=field;
end

V=nc_vnames(ncfile);
nvars=length(V.Variables);

for n=1:nvars
  name=char(V.Variables(n).Name);
  switch name
    case {Vname.mask_rho}
      got.mask_rho=true;
      Vid.(name)=netcdf.inqVarID(ncid,name);
    case {Vname.mask_psi}
      got.mask_psi=true;
      Vid.(name)=netcdf.inqVarID(ncid,name);
   case {Vname.mask_u}
      got.mask_u  =true;
      Vid.(name)=netcdf.inqVarID(ncid,name);
    case {Vname.mask_v}
      got.mask_v  =true;
      Vid.(name)=netcdf.inqVarID(ncid,name);
  end
end

%--------------------------------------------------------------------------
% If appropriate, open GRID NetCDF and define Land/Sea mask varaibles.
%--------------------------------------------------------------------------

defmode=~got.mask_rho || ~got.mask_psi || ~got.mask_u || ~got.mask_v;

if (defmode)

% Inquire about dimensions.

  Dname.xr='xi_rho';  Dname.yr='eta_rho';
  Dname.xp='xi_psi';  Dname.yp='eta_psi';
  Dname.xu='xi_u';    Dname.yu='eta_u';
  Dname.xv='xi_v';    Dname.yv='eta_v';

  D=nc_dinfo(ncfile);
  ndims=length(D);
  for n=1:ndims
    name=char(D(n).Name);
    switch name
      case {Dname.xr}
        Dsize.xr=D(n).Length;
        did.xr  =D(n).dimid;
      case {Dname.yr}
        Dsize.yr=D(n).Length;
        did.yr  =D(n).dimid;
      case {Dname.xp}
        Dsize.xp=D(n).Length;
        did.xp  =D(n).dimid;
      case {Dname.yp}
        Dsize.yp=D(n).Length;
        did.yp  =D(n).dimid;
      case {Dname.xu}
        Dsize.xu=D(n).Length;
        did.xu  =D(n).dimid;
      case {Dname.yu}
        Dsize.yu=D(n).Length;
        did.yu  =D(n).dimid;
      case {Dname.xv}
        Dsize.xv=D(n).Length;
        did.xv  =D(n).dimid;
      case {Dname.yv}
        Dsize.yv=D(n).Length;
        did.yv  =D(n).dimid;
    end
  end

% Put GRID NetCDF file in definition mode.

  netcdf.reDef(ncid);

% Define RHO-points mask.

  if (~got.mask_rho)
    field='mask_rho';
    Vid.(field)=netcdf.defVar(ncid, field, 'NC_DOUBLE', [did.xr did.yr]);
    netcdf.putAtt(ncid, Vid.(field), 'long_name',                       ...
                  'land/sea mask on RHO-points');
    netcdf.putAtt(ncid, Vid.(field), 'flag_values',                     ...
                  [double(0) double(1)]);
    netcdf.putAtt(ncid, Vid.(field), 'flag_meanings',                   ...
                  'land water');
    netcdf.putAtt(ncid, Vid.(field), 'coordinates',                     ...
                  'lon_rho lat_rho');
  end

% Define PSI-points mask.

  if (~got.mask_psi)
    field='mask_psi';
    Vid.(field)=netcdf.defVar(ncid, field, 'NC_DOUBLE', [did.xp did.yp]);
    netcdf.putAtt(ncid, Vid.(field), 'long_name',                       ...
                  'land/sea mask on PSI-points');
    netcdf.putAtt(ncid, Vid.(field), 'flag_values',                     ...
                  [double(0) double(1)]);
    netcdf.putAtt(ncid, Vid.(field), 'flag_meanings',                   ...
                  'land water');
    netcdf.putAtt(ncid, Vid.(field), 'coordinates',                     ...
                  'lon_psi lat_psi');
  end

% Define U-points mask.

  if (~got.mask_u)
    field='mask_u';
    Vid.(field)=netcdf.defVar(ncid, field, 'NC_DOUBLE', [did.xu did.yu]);
    netcdf.putAtt(ncid, Vid.(field), 'long_name',                       ...
                  'land/sea mask on U-points');
    netcdf.putAtt(ncid, Vid.(field), 'flag_values',                     ...
                  [double(0) double(1)]);
    netcdf.putAtt(ncid, Vid.(field), 'flag_meanings',                   ...
                  'land water');
    netcdf.putAtt(ncid, Vid.(field), 'coordinates',                     ...
                  'lon_u lat_u');
  end

% Define V-points mask.

  if (~got.mask_v)
    field='mask_v';
    Vid.(field)=netcdf.defVar(ncid, field, 'NC_DOUBLE', [did.xv did.yv]);
    netcdf.putAtt(ncid, Vid.(field), 'long_name',                       ...
                  'land/sea mask on V-points');
    netcdf.putAtt(ncid, Vid.(field), 'flag_values',                     ...
                  [double(0) double(1)]);
    netcdf.putAtt(ncid, Vid.(field), 'flag_meanings',                   ...
                  'land water');
    netcdf.putAtt(ncid, Vid.(field), 'coordinates',                     ...
                  'lon_v lat_v');
  end

% Leave definition mode.

  netcdf.endDef(ncid)

end

%--------------------------------------------------------------------------
% Get grid coordinates at rho-points.
%--------------------------------------------------------------------------

Vid.lon_rho = netcdf.inqVarID(ncid, 'lon_rho');
F.lon_rho = netcdf.getVar(ncid, Vid.lon_rho, 'double');

Vid.lat_rho = netcdf.inqVarID(ncid, 'lat_rho');
F.lat_rho = netcdf.getVar(ncid, Vid.lat_rho, 'double');

lon_rho=F.lon_rho';
lat_rho=F.lat_rho';

%--------------------------------------------------------------------------
% Extract coastlines from GSHHS database.
%--------------------------------------------------------------------------

if (extract_coast)

  dx=5*abs(mean(mean(gradient(lon_rho))));
  dy=5*abs(mean(mean(gradient(lat_rho))));

  F.Llon=min(min(lon_rho));  F.Llon=F.Llon-dx;
  F.Rlon=max(max(lon_rho));  F.Rlon=F.Rlon+dx;
  F.Blat=min(min(lat_rho));  F.Blat=F.Blat-dy;
  F.Tlat=max(max(lat_rho));  F.Tlat=F.Tlat+dy;

  [C]=r_gshhs(F.Llon,F.Rlon,F.Blat,F.Tlat,Cname);
  [C]=x_gshhs(F.Llon,F.Rlon,F.Blat,F.Tlat,C,'patch');
  clon=C.lon';
  clat=C.lat';

  F.clon=clon;
  F.clat=clat;

  clear C

end

%--------------------------------------------------------------------------
% Compute Land/Sea mask on RHO-points.
%--------------------------------------------------------------------------

[i,~]=size(clon);  if (i > 1), clon=clon'; end
[i,~]=size(clat);  if (i > 1), clat=clat'; end

f=find(~isfinite(clon) | ~isfinite(clat));
f=f(:).';
if ~any(f)
  f=[0 length(clon)+1];
end
if (f(1) ~= 1)
  f=[0 f];
end
if (f(end) ~= length(clon))
  f(end+1)=length(clon)+1;
end

Mask=zeros(size(lon_rho));

for i=2:length(f)
  g=find(Mask == 0);
  if (~any(g))
    break
  end
  j=f(i-1)+1:f(i)-1;
  if (length(j) > 2)
    Mask(g)=inpolygon(lon_rho(g),lat_rho(g),clon(j),clat(j));
  end
end

mask_rho=1-Mask;
F.mask_rho=mask_rho';

clear Mask mask_rho clon clat f g lon_rho lat_rho

%--------------------------------------------------------------------------
% Compute Land/Sea mask on PSI-, U-, and V-points.
%--------------------------------------------------------------------------

[F.mask_u, F.mask_v, F.mask_psi]=uvp_masks(F.mask_rho);

%--------------------------------------------------------------------------
% Write out Land/Sea masking into GRID NetCDF file.
%--------------------------------------------------------------------------

disp(' ');
for var = VarList
  field=char(var);
  disp(['writing ',sprintf('%10s',field),                               ...
        ',   size = ', num2str(size(F.(field)))]);
  netcdf.putVar(ncid, Vid.(field), F.(field));
end

netcdf.close(ncid)

disp(' ');
F.profile=strcat(num2str(toc(Tstart)),' sec');
disp(['Elapsed time to compute mask ', F.profile])

return
