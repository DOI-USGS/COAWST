function [F]=landsea(Gname, clon, clat);

% LANDSEA:  Computes and writes ROMS Land/Sea masks
%
% [F]=landsea(Gname,Clon,Clat)
%
% This function computes the Land/Sea mask for an application grid and
% writes it into input NetCDF file. The Land/Sea mask is computed from
% closed coastlines polygoms extracted from the GSHHS dataset.
%
% On Input:
%
%    Gname       NetCDF file name (character string).
%    Clon        Coastline longitude (real vector).
%    Clat        Coastline latitude  (real vector).
%
% On Output:
%
%    F           Land/Sea mask (structure array):
%                  F.rmask =>  Mask at rho-points.
%                  F.pmask =>  Mask at psi-points.
%                  F.umask =>  Mask at u-points.
%                  F.vmask =>  Mask at v-points.
%
% Calls:   MEXNC (Interface to NetCDF library using Matlab):
%          nc_dim, nc_read, nc_vdef, nc_vname, nc_write
%
%          r_gshhs, x_gshhs, mexinside
%
% The masking algorithm is adapted from Charles R. Denham routine
% "domask.m" from the SeaGrid package.
%

% svn $Id: landsea.m 485 2010-07-07 18:10:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

database=0;
extract_coast=0;
DIR='C:\Documents and Settings\cdenamie\My Documents\MATLAB\ROMS_tools\GSHHS\';  % GSHHS directory path

if (nargin < 3),
  extract_coast=1;              % Extract coastlines from GSHHS database
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
  end,

  Cname=strcat(DIR,name);
  F.Cname=Cname;

end,

F.Gname=Gname;

%-----------------------------------------------------------------------
% Inquire grid NetCDF file about mask variables.
%-----------------------------------------------------------------------

got.rmask=0;  Vname.rmask='mask_rho';
got.pmask=0;  Vname.pmask='mask_psi';
got.umask=0;  Vname.umask='mask_u';
got.vmask=0;  Vname.vmask='mask_v';

[varnam,nvars]=nc_vname(Gname);
for n=1:nvars,
  name=deblank(varnam(n,:));
  switch name
    case {Vname.rmask}
      got.rmask=1;
    case {Vname.pmask}
      got.pmask=1;
    case {Vname.umask}
      got.umask=1;
    case {Vname.vmask}
      got.vmask=1;
  end,
end,

%-----------------------------------------------------------------------
% If appropriate, open GRID NetCDF and define Land/Sea mask varaibles.
%-----------------------------------------------------------------------

defmode=~got.rmask | ~got.pmask | ~got.umask | ~got.vmask;

if (defmode),

% Inquire about dimensions.

  Dname.xr='xi_rho';  Dname.yr='eta_rho';
  Dname.xp='xi_psi';  Dname.yp='eta_psi';
  Dname.xu='xi_u';    Dname.yu='eta_u';
  Dname.xv='xi_v';    Dname.yv='eta_v';

  [Dnames,Dsizes]=nc_dim(Fname);
  ndims=length(Dsizes);
  for n=1:ndims,
    dimid=n;
    name=deblank(Dnames(n,:));
    switch name
      case {Dname.xr}
        Dsize.xr=Dsizes(n);
        did.xr=dimid-1;
      case {Dname.yr}
        Dsize.yr=Dsizes(n);
        did.yr=dimid-1;
      case {Dname.xp}
        Dsize.xp=Dsizes(n);
        did.xp=dimid-1;
      case {Dname.yp}
        Dsize.yp=Dsizes(n);
        did.yp=dimid-1;
      case {Dname.xu}
        Dsize.xu=Dsizes(n);
        did.xu=dimid-1;
      case {Dname.yu}
        Dsize.yu=Dsizes(n);
        did.yu=dimid-1;
      case {Dname.xv}
        Dsize.xv=Dsizes(n);
        did.xv=dimid-1;
      case {Dname.yv}
        Dsize.yv=Dsizes(n);
        did.yv=dimid-1;
    end,
  end,

% Open GRID NetCDF file.

  [ncid,status]=mexnc('open',Gname,'nc_write');
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['LANDSEA: OPEN - unable to open file: ', Gname])
    return
  end

% Put GRID NetCDF file in definition mode.

  [ncglobal]=mexnc('parameter','nc_global');
  [ncdouble]=mexnc('parameter','nc_double');
  [ncunlim ]=mexnc('parameter','nc_unlimited');
  [ncfloat ]=mexnc('parameter','nc_float');
  [ncchar  ]=mexnc('parameter','nc_char');

  [status]=mexnc('redef',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['LANDSEA: REDEF - unable to put into define mode.'])
    return
  end

% Define RHO-points mask.

  if (~got.rmask),
    Var.name          = Vname.rmask;
    Var.type          = ncfloat;
    Var.dimid         = [did.yr did.xr];
    Var.long_name     = 'mask on RHO-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1), ...
                         'water'];
    [varid,status]=nc_vdef(ncid,Var);
    clear Var
  end,

% Define PSI-points mask.

  if (~got.pmask),
    Var.name          = Vname.pmask;
    Var.type          = ncfloat;
    Var.dimid         = [did.yp did.xp];
    Var.long_name     = 'mask on PSI-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1), ...
                         'water'];
    [varid,status]=nc_vdef(ncid,Var);
    clear Var
  end,

% Define U-points mask.

  if (~got.umask),
    Var.name          = Vname.umask;
    Var.type          = ncfloat;
    Var.dimid         = [did.yu did.xu];
    Var.long_name     = 'mask on U-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1), ...
                         'water'];
    [varid,status]=nc_vdef(ncid,Var);
    clear Var
  end,

% Define V-points mask.

  if (~got.vmask),
    Var.name          = Vname.vmask;
    Var.type          = ncfloat;
    Var.dimid         = [did.yv did.xv];
    Var.long_name     = 'mask on V-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1), ...
                         'water'];
    [varid,status]=nc_vdef(ncid,Var);
    clear Var
  end,

% Leave definition mode.

  [status]=mexnc('enddef',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['LANDSEA: ENDDEF - unable to leave definition mode.'])
  end,

% Close GRID NetCDF file.

  [status]=mexnc('close',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['LANDSEA: CLOSE - unable to close NetCDF file: ', Gname])
  end,

end,

%-----------------------------------------------------------------------
% Get grid coordinates at rho-points.
%-----------------------------------------------------------------------

F.rlon=nc_read(Gname,'lon_rho');
F.rlat=nc_read(Gname,'lat_rho');

rlon=F.rlon';
rlat=F.rlat';
[Jm,Im]=size(rlon);

%-----------------------------------------------------------------------
% Extract coastlines from GSHHS database.
%-----------------------------------------------------------------------

if (extract_coast),

  dx=5*abs(mean(mean(gradient(rlon))));
  dy=5*abs(mean(mean(gradient(rlat))));

  F.Llon=min(min(rlon));  F.Llon=F.Llon-dx;
  F.Rlon=max(max(rlon));  F.Rlon=F.Rlon+dx;
  F.Blat=min(min(rlat));  F.Blat=F.Blat-dy;
  F.Tlat=max(max(rlat));  F.Tlat=F.Tlat+dy;

  h=warndlg('Please wait ...', 'Extracting Coastlines');
  drawnow

  [C]=r_gshhs(F.Llon,F.Rlon,F.Blat,F.Tlat,Cname);
  [C]=x_gshhs(F.Llon,F.Rlon,F.Blat,F.Tlat,C,'patch');
  clon=C.lon';
  clat=C.lat';

  F.clon=clon;
  F.clat=clat;

  if (ishandle(h)),
    delete(h),
  end,

  clear C

end,

%-----------------------------------------------------------------------
% Compute Land/Sea mask on RHO-points.
%-----------------------------------------------------------------------

[i,j]=size(clon);  if (i > 1), clon=clon'; end,
[i,j]=size(clat);  if (i > 1), clat=clat'; end,

f=find(~isfinite(clon) | ~isfinite(clat));
f=f(:).';
if ~any(f),
  f=[0 length(clon)+1];
end,
if (f(1) ~= 1),
  f=[0 f];
end,
if (f(end) ~= length(clon)),
  f(end+1)=length(clon)+1;
end

h=warndlg('Please wait ...', 'Computing Mask');
drawnow

Mask=zeros(size(rlon));

for i=2:length(f),
  g=find(Mask == 0);
  if (~any(g)),
    break,
  end,
  j=f(i-1)+1:f(i)-1;
  if (length(j) > 2),
    Mask(g)=mexinside(rlon(g),rlat(g),clon(j),clat(j));
  end,
end,

if (ishandle(h)),
  delete(h),
end,

rmask=1-Mask;

clear Mask clon clat f g rlon rlat

%-----------------------------------------------------------------------
% Compute Land/Sea mask on PSI-, U-, and V-points.
%-----------------------------------------------------------------------

pmask(1:Jm-1,1:Im-1)=rmask(2:Jm  ,2:Im  ).* ...
                     rmask(2:Jm  ,1:Im-1).* ...
                     rmask(1:Jm-1,2:Im  ).* ...
                     rmask(1:Jm-1,1:Im-1);
umask(1:Jm  ,1:Im-1)=rmask(1:Jm  ,2:Im  ).* ...
                     rmask(1:Jm  ,1:Im-1);
vmask(1:Jm-1,1:Im  )=rmask(2:Jm  ,1:Im  ).* ...
                     rmask(1:Jm-1,1:Im  );

%-----------------------------------------------------------------------
% Write out Land/Sea masking into GRID NetCDF file.
%-----------------------------------------------------------------------

F.rmask=rmask';  clear rmask
F.pmask=pmask';  clear pmask
F.umask=umask';  clear umask
F.vmask=vmask';  clear vmask

[status]=nc_write(Gname,Vname.rmask,F.rmask);
[status]=nc_write(Gname,Vname.pmask,F.pmask);
[status]=nc_write(Gname,Vname.umask,F.umask);
[status]=nc_write(Gname,Vname.vmask,F.vmask);

return


