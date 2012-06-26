function status=write_bath(Gname, rbath);

% WRITE_SCOPE:  Writes ROMS bathymetry
%
% status=write_mask(Gname,rbath)
%
% This routine writes out mask data into GRID NetCDF file.
%
% On Input:
%
%    Gname       GRID NetCDF file name (character string).
%    rbath       Bathymetry on RHO-points (real matrix):
%

% svn $Id: write_mask.m 485 2010-07-07 18:10:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%---------------------------------------------------------------------------
% Inquire grid NetCDF file about mask variables.
%--------------------------------------------------------------------------
got.rbath=0;  define.rbath=1;  Vname.rbath='h';

[varnam,nvars]=nc_vname(Gname);
for n=1:nvars
  name=deblank(varnam(n,:));
  switch name
    case {Vname.rbath}
      got.rbath=1;
      define.rbath=0;
  end
end
%---------------------------------------------------------------------------
%  If appropriate, define Land/Sea mask variables.
%--------------------------------------------------------------------------
if (define.rbath)
%  Inquire about dimensions.
  Dname.xp='xi_psi';    Dname.yp='eta_psi';
  Dname.xr='xi_rho';    Dname.yr='eta_rho';
  Dname.xu='xi_u';      Dname.yu='eta_u';
  Dname.xv='xi_v';      Dname.yv='eta_v';
  [Dnames,Dsizes]=nc_dim(Gname);
  ndims=length(Dsizes);
  for n=1:ndims
    dimid=n-1;
    name=deblank(Dnames(n,:));
    switch name
      case {Dname.xp}
        Dsize.xp=Dsizes(n);
        did.xp=dimid;
      case {Dname.yp}
        Dsize.yp=Dsizes(n);
        did.yp=dimid;
      case {Dname.xr}
        Dsize.xr=Dsizes(n);
        did.xr=dimid;
      case {Dname.yr}
        Dsize.yr=Dsizes(n);
        did.yr=dimid;
      case {Dname.xu}
        Dsize.xu=Dsizes(n);
        did.xu=dimid;
      case {Dname.yu}
        Dsize.yu=Dsizes(n);
        did.yu=dimid;
      case {Dname.xv}
        Dsize.xv=Dsizes(n);
        did.xv=dimid;
      case {Dname.yv}
        Dsize.yv=Dsizes(n);
        did.yv=dimid;
    end
  end
%  Define NetCDF parameters.
  [ncglobal]=mexnc('parameter','nc_global');
  [ncdouble]=mexnc('parameter','nc_double');
  [ncfloat ]=mexnc('parameter','nc_float');
  [ncchar  ]=mexnc('parameter','nc_char');
%  Open GRID NetCDF file.
  [ncid,status]=mexnc('open',Gname,'nc_write');
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['WRITE_BATH: OPEN - unable to open file: ', Gname]);
  end

%  Put GRID NetCDF file into define mode.
  [status]=mexnc('redef',ncid);
  if (status ~= 0)
    disp('  ');
    disp(mexnc('strerror',status));
    error('WRITE_BATH: REDEF - unable to put into define mode.');
  end

%  Define Bathymetry on RHO-points.
  if (define.rbath)
    Var.name        = Vname.rbath;
    Var.type        = ncdouble;
    Var.dimid       = [did.yr did.xr];
    Var.long_name   = 'Final bathymetry at RHO-points';
    Var.units       = 'meter';
    Var.coordinates = 'lon_rho lat_rho';
    Var.field       = 'bath, scalar';
    [~,status]      = nc_vdef(ncid,Var);
    clear Var
  end

%  Leave definition mode.
  [status]=mexnc('enddef',ncid);
  if (status ~= 0)
    disp('  ');
    disp(mexnc('strerror',status));
    error('WRITE_BATH: ENDDEF - unable to leave definition mode.');
  end
%  Close GRID NetCDF file.
  [status]=mexnc('close',ncid);
  if (status ~= 0)
    disp('  ');
    disp(mexnc('strerror',status));
    error(['WRITE_BATH: CLOSE - unable to close NetCDF file: ', Gname]);
  end
end

%--------------------------------------------------------------------------
%  Write out bathymetry data into GRID NetCDF file.
%--------------------------------------------------------------------------
[status]=nc_write(Gname,Vname.rbath,rbath);
return
