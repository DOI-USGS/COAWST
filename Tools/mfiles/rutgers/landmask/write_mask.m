function status = write_mask(ncfile, rmask, umask, vmask, pmask)

% WRITE_MASK:  Writes ROMS Land/Sea masks
%
% status = write_mask(ncfile, rmask, umask, vmask, pmask)
%
% This routine writes out mask data into GRID NetCDF file.
%
% On Input:
%
%    ncfile      GRID NetCDF file name (character string).
%    rmask       Land/Sea mask on RHO-points (real matrix):
%                  rmask=0 land, rmask=1 Sea
%    umask       Land/Sea mask on U-points (real matrix):
%                  umask=0 land, umask=1 Sea
%    vmask       Land/Sea mask on V-points (real matrix):
%                  vmask=0 land, vmask=1 Sea
%    pmask       Land/Sea mask on PSI-points (real matrix):
%                  pmask=0 land, pmask=1 Sea
%

% svn $Id: write_mask.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire grid NetCDF file about mask variables.
%--------------------------------------------------------------------------

got.pmask=false;  define.pmask=true;  Vname.pmask='mask_psi';
got.rmask=false;  define.rmask=true;  Vname.rmask='mask_rho';
got.umask=false;  define.umask=true;  Vname.umask='mask_u';
got.vmask=false;  define.vmask=true;  Vname.vmask='mask_v';

V=nc_vnames(ncfile);
nvars=length(V.Variables);
for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch name
    case {Vname.pmask}
      got.pmask   =true;
      define.pmask=false;
    case {Vname.rmask}
      got.rmask   =true;
      define.rmask=false;
    case {Vname.umask}
      got.umask   =true;
      define.umask=false;
    case {Vname.vmask}
      got.vmask   =true;
      define.vmask=false;
  end
end

%-------------------------------------------------------------------------
%  If appropriate, define Land/Sea mask variables.
%-------------------------------------------------------------------------

if (define.pmask || define.rmask || define.umask || define.pmask),

%  Inquire about dimensions.

  Dname.xp='xi_psi';    Dname.yp='eta_psi';
  Dname.xr='xi_rho';    Dname.yr='eta_rho';
  Dname.xu='xi_u';      Dname.yu='eta_u';
  Dname.xv='xi_v';      Dname.yv='eta_v';

  D=nc_dinfo(ncfile);
  ndims=length(D);
  for n=1:ndims,
    name=char(D(n).Name);
    switch name
      case {Dname.xp}
        Dsize.xp=D(n).Length;
        did.xp  =D(n).dimid;
      case {Dname.yp}
        Dsize.yp=D(n).Length;
        did.yp  =D(n).dimid;
      case {Dname.xr}
        Dsize.xr=D(n).Length;
        did.xr  =D(n).dimid;
      case {Dname.yr}
        Dsize.yr=D(n).Length;
        did.yr  =D(n).dimid;
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
    end,
  end,

%  Open GRID NetCDF file.

  [ncid,status]=mexnc('open',ncfile,'nc_write');
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['WRITE_MASK: OPEN - unable to open file: ', ncfile]);
  end


%  Put GRID NetCDF file into define mode.

  [status]=mexnc('redef',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error('WRITE_MASK: REDEF - unable to put into define mode.');
  end

%  Define Land/Sea mask on RHO-points.

  if (define.rmask),
    Var.name          = Vname.rmask;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yr did.xr];
    Var.long_name     = 'mask on RHO-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1), ...
                         'water'];

    [~,status]=nc_vdef(ncid,Var);
    clear Var
  end

%  Define Land/Sea mask on PSI-points.

  if (define.pmask),
    Var.name          = Vname.pmask;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yp did.xp];
    Var.long_name     = 'mask on PSI-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1), ...
                         'water'];

    [~,status]=nc_vdef(ncid,Var);
    clear Var
  end

%  Define Land/Sea mask on U-points.

  if (define.umask),
    Var.name          = Vname.umask;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yu did.xu];
    Var.long_name     = 'mask on U-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1), ...
                         'water'];

    [~,status]=nc_vdef(ncid,Var);
    clear Var
  end

%  Define Land/Sea mask on V-points.

  if (define.vmask),
    Var.name          = Vname.vmask;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yv did.xv];
    Var.long_name     = 'mask on V-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['land', blanks(1), ...
                         'water'];

    [~,status]=nc_vdef(ncid,Var);
    clear Var
  end

%  Leave definition mode.

  [status]=mexnc('enddef',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error('WRITE_MASK: ENDDEF - unable to leave definition mode.');
  end

%  Close GRID NetCDF file.

  [status]=mexnc('close',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['WRITE_MASK: CLOSE - unable to close NetCDF file: ', ncfile]);
  end

end

%--------------------------------------------------------------------------
%  Write out mask data into GRID NetCDF file.
%--------------------------------------------------------------------------

[status]=nc_write(ncfile,Vname.rmask,rmask);
[status]=nc_write(ncfile,Vname.pmask,pmask);
[status]=nc_write(ncfile,Vname.umask,umask);
[status]=nc_write(ncfile,Vname.vmask,vmask);

return
