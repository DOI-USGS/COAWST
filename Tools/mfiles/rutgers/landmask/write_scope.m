function status = write_scope(ncfile, Rscope, Uscope, Vscope)

% WRITE_SCOPE:  Writes ROMS adjoint sensitivity scope masks
%
% status = write_scope(ncfile, Rscope, Uscope, Vscope)
%
% This routine writes out adjoint sensitivity scope mask data into
% GRID NetCDF file.
%
% On Input:
%
%    ncfile      GRID NetCDF file name (character string)
%    Rscope      Adjoint sensitivity scope mask on RHO-points:
%                  Rscope=0 inactive, Rscope=1 active
%    Uscope      Adjoint sensitivity scope mask on U-points:
%                  Uscope=0 inactive, Uscope=1 active
%    Vscope      Adjoint sensitivity scope mask on V-points:
%                  Vscope=0 inactive, Vscope=1 active
%

% svn $Id: write_scope.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%-------------------------------------------------------------------------
% Inquire grid NetCDF file about scope variables.
%-------------------------------------------------------------------------

got.Rscope=false;  define.Rscope=true;  Vname.Rscope='scope_rho';
got.Uscope=false;  define.Uscope=true;  Vname.Uscope='scope_u';
got.Vscope=false;  define.Vscope=true;  Vname.Vscope='scope_v';

V=nc_vnames(ncfile);
nvars=length(V.Variables);
for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch name
    case {Vname.Rscope}
      got.Rscope   =true;
      define.Rscope=false;
    case {Vname.Uscope}
      got.Uscope   =true;
      define.Uscope=false;
    case {Vname.Vscope}
      got.Vscope   =true;
      define.Vscope=false;
  end,
end,

%--------------------------------------------------------------------------
%  If appropriate, define scope variables.
%--------------------------------------------------------------------------

if (define.Rscope || define.Uscope || define.Vscope),

%  Inquire about dimensions.

  Dname.xr='xi_rho';    Dname.yr='eta_rho';
  Dname.xu='xi_u';      Dname.yu='eta_u';
  Dname.xv='xi_v';      Dname.yv='eta_v';

  D=nc_dinfo(ncfile);
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
    error('WRITE_MASK: REFDEF - unable to put into define mode.');
  end

%  Define scope mask on RHO-points.

  if (define.Rscope),
    Var.name          = Vname.Rscope;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yr did.xr];
    Var.long_name     = 'adjoint sensitivity scope mask on RHO-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['inactive', blanks(1), ...
                         'active'];

    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end

%  Define scope mask on U-points.

  if (define.Uscope),
    Var.name          = Vname.Uscope;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yu did.xu];
    Var.long_name     = 'adjoint sensitivity scope mask on U-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['inactive', blanks(1), ...
                         'active'];

    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end

%  Define scope mask on V-points.

  if (define.Vscope),
    Var.name          = Vname.Vscope;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.yv did.xv];
    Var.long_name     = 'adjoint sensitivity scope mask on V-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['inactive', blanks(1), ...
                         'active'];

    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
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

[status]=nc_write(ncfile,Vname.Rscope,Rscope);
[status]=nc_write(ncfile,Vname.Uscope,Uscope);
[status]=nc_write(ncfile,Vname.Vscope,Vscope);

return
