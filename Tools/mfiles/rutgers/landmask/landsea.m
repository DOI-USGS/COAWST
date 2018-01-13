function [F]=landsea(ncfile, clon, clat)

% LANDSEA:  Computes and writes ROMS Land/Sea masks
%
% [F]=landsea(ncfile,Clon,Clat)
%
% Computes automatically the Land/Sea mask for an application grid and
% writes it into input NetCDF file. The Land/Sea mask is computed from
% closed coastlines polygoms extracted from the GSHHS dataset.
%
% On Input:
%
%    ncfile      NetCDF file name (string)
%    Clon        Coastline longitude (real vector)
%    Clat        Coastline latitude  (real vector)
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

% svn $Id: landsea.m 832 2017-01-24 22:07:36Z arango $
%=========================================================================%
%  Copyright (c) 2002-2017 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

extract_coast=false;

if (nargin < 3),
  DIR='~/ocean/GSHHS/Version_1.2';                % GSHHS directory path

  extract_coast=true;           % Extract coastlines from GSHHS database
% database='full';              % Full resolution database
% database='high';              % High resolution database
  database='intermediate';      % Intermediate resolution database
% database='low';               % Low resolution database
% database='crude';             % crude resolution database

  switch database,
    case 'full'
      name='gshhs_f.b';
    case 'high'
      name='gshhs_h.b';
    case 'intermediate'
      name='gshhs_i.b';
    case 'low'
      name='gshhs_l.b';
    case 'crude'
      name='gshhs_c.b';
  end

  Cname=fullfile(DIR,name);
  F.Cname=Cname;
end

F.ncfile=ncfile;

%-----------------------------------------------------------------------
% Inquire grid NetCDF file about mask variables.
%-----------------------------------------------------------------------

got.mask_rho=false;  Vname.mask_rho='mask_rho';
got.mask_psi=false;  Vname.mask_psi='mask_psi';
got.mask_u  =false;  Vname.mask_u  ='mask_u';
got.mask_v  =false;  Vname.mask_v  ='mask_v';

V=nc_vnames(ncfile);
nvars=length(V.Variables);
for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch name
    case {Vname.mask_rho}
      got.mask_rho=true;
    case {Vname.mask_psi}
      got.mask_psi=true;
    case {Vname.mask_u}
      got.mask_u  =true;
    case {Vname.mask_v}
      got.mask_v  =true;
  end
end

%--------------------------------------------------------------------------
% If appropriate, open GRID NetCDF and define Land/Sea mask varaibles.
%--------------------------------------------------------------------------

defmode=~got.mask_rho || ~got.mask_psi || ~got.mask_u || ~got.mask_v;

if (defmode),

% Inquire about dimensions.

  Dname.xr='xi_rho';  Dname.yr='eta_rho';
  Dname.xp='xi_psi';  Dname.yp='eta_psi';
  Dname.xu='xi_u';    Dname.yu='eta_u';
  Dname.xv='xi_v';    Dname.yv='eta_v';

  D=nc_dinfo(Fname);
  ndims=length(D);
  for n=1:ndims,
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

% Open GRID NetCDF file.

  [ncid,status]=mexnc('open',ncfile,'nc_write');
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['LANDSEA: OPEN - unable to open file: ', ncfile])
  end

% Put GRID NetCDF file in definition mode.

  [status]=mexnc('redef',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error('LANDSEA: REDEF - unable to put into define mode.')
  end

% Define RHO-points mask.

  if (~got.mask_rho),
    Var.name          = Vname.mask_rho;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yr did.xr];
    Var.long_name     = 'mask on RHO-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1),                             ...
                         'water'];
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end

% Define PSI-points mask.

  if (~got.mask_psi),
    Var.name          = Vname.mask_psi;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yp did.xp];
    Var.long_name     = 'mask on PSI-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1),                             ...
                         'water'];
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end

% Define U-points mask.

  if (~got.mask_u),
    Var.name          = Vname.mask_u;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yu did.xu];
    Var.long_name     = 'mask on U-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1),                             ...
                         'water'];
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end

% Define V-points mask.

  if (~got.mask_v),
    Var.name          = Vname.mask_v;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yv did.xv];
    Var.long_name     = 'mask on V-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1), ...
                         'water'];
    [~,status]=nc_vdef(ncid,Var);
    clear Var
    if (status ~= 0), return, end
  end

% Leave definition mode.

  [status]=mexnc('enddef',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error('LANDSEA: ENDDEF - unable to leave definition mode.')
  end

% Close GRID NetCDF file.

  [status]=mexnc('close',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['LANDSEA: CLOSE - unable to close NetCDF file: ', ncfile])
  end

end

%--------------------------------------------------------------------------
% Get grid coordinates at rho-points.
%--------------------------------------------------------------------------

F.lon_rho=nc_read(ncfile,'lon_rho');
F.lat_rho=nc_read(ncfile,'lat_rho');

lon_rho=F.lon_rho';
lat_rho=F.lat_rho';
[Jm,Im]=size(lon_rho);

%--------------------------------------------------------------------------
% Extract coastlines from GSHHS database.
%--------------------------------------------------------------------------

if (extract_coast),

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

[i,~]=size(clon);  if (i > 1), clon=clon'; end,
[i,~]=size(clat);  if (i > 1), clat=clat'; end,

f=find(~isfinite(clon) | ~isfinite(clat));
f=f(:).';
if ~any(f),
  f=[0 length(clon)+1];
end
if (f(1) ~= 1),
  f=[0 f];
end
if (f(end) ~= length(clon)),
  f(end+1)=length(clon)+1;
end

Mask=zeros(size(lon_rho));

for i=2:length(f),
  g=find(Mask == 0);
  if (~any(g)),
    break
  end
  j=f(i-1)+1:f(i)-1;
  if (length(j) > 2),
    Mask(g)=inpolygon(lon_rho(g),lat_rho(g),clon(j),clat(j));
  end
end

mask_rho=1-Mask;

clear Mask clon clat f g lon_rho lat_rho

%--------------------------------------------------------------------------
% Compute Land/Sea mask on PSI-, U-, and V-points.
%--------------------------------------------------------------------------

mask_psi(1:Jm-1,1:Im-1)=mask_rho(2:Jm  ,2:Im  ).*                       ...
                        mask_rho(2:Jm  ,1:Im-1).*                       ...
                        mask_rho(1:Jm-1,2:Im  ).*                       ...
                        mask_rho(1:Jm-1,1:Im-1);
mask_u  (1:Jm  ,1:Im-1)=mask_rho(1:Jm  ,2:Im  ).*                       ...
                        mask_rho(1:Jm  ,1:Im-1);
mask_v  (1:Jm-1,1:Im  )=mask_rho(2:Jm  ,1:Im  ).*                       ...
                        mask_rho(1:Jm-1,1:Im  );

%--------------------------------------------------------------------------
% Write out Land/Sea masking into GRID NetCDF file.
%--------------------------------------------------------------------------

F.mask_rho=mask_rho';
F.mask_psi=mask_psi';
F.mask_u  =mask_u';
F.mask_v  =mask_v';

[status]=nc_write(ncfile,Vname.mask_rho,F.mask_rho);
if (status ~= 0), return, end

[status]=nc_write(ncfile,Vname.mask_psi,F.mask_psi);
if (status ~= 0), return, end

[status]=nc_write(ncfile,Vname.mask_u  ,F.mask_u  );
if (status ~= 0), return, end

[status]=nc_write(ncfile,Vname.mask_v  ,F.mask_v  );
if (status ~= 0), return, end

return


