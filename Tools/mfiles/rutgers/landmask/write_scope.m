function status=write_scope(Gname, Rscope, Uscope, Vscope);

% WRITE_SCOPE:  Writes ROMS adjoint sensitivity scope masks
%
% status=write_scope(Gname, Rscope, Uscope, Vscope)
%
% This routine writes out adjoint sensitivity scope mask data into
% GRID NetCDF file.
%
% On Input:
%
%    Gname       GRID NetCDF file name (character string).
%    Rscope      Adjoint sensitivity scope mask on RHO-points:
%                  Rscope=0 inactive, Rscope=1 active.
%    Uscope      Adjoint sensitivity scope mask on U-points:
%                  Uscope=0 inactive, Uscope=1 active.
%    Vscope      Adjoint sensitivity scope mask on V-points:
%                  Vscope=0 inactive, Vscope=1 active.
%

% svn $Id: write_scope.m 485 2010-07-07 18:10:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%---------------------------------------------------------------------------
% Inquire grid NetCDF file about scope variables.
%---------------------------------------------------------------------------

got.Rscope=0;  define.Rscope=1;  Vname.Rscope='scope_rho';
got.Uscope=0;  define.Uscope=1;  Vname.Uscope='scope_u';
got.Vscope=0;  define.Vscope=1;  Vname.Vscope='scope_v';

[varnam,nvars]=nc_vname(Gname);
for n=1:nvars,
  name=deblank(varnam(n,:));
  switch name
    case {Vname.Rscope}
      got.Rscope=1;
      define.Rscope=0;
    case {Vname.Uscope}
      got.Uscope=1;
      define.Uscope=0;
    case {Vname.Vscope}
      got.Vscope=1;
      define.Vscope=0;
  end,
end,

%---------------------------------------------------------------------------
%  If appropriate, define scope variables.
%---------------------------------------------------------------------------

if (define.Rscope | define.Uscope | define.Vscope),

%  Inquire about dimensions.

  Dname.xr='xi_rho';    Dname.yr='eta_rho';
  Dname.xu='xi_u';      Dname.yu='eta_u';
  Dname.xv='xi_v';      Dname.yv='eta_v';

  [Dnames,Dsizes]=nc_dim(Gname);
  ndims=length(Dsizes);
  for n=1:ndims,
    dimid=n-1;
    name=deblank(Dnames(n,:));
    switch name
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
    end,
  end,

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
    error(['WRITE_MASK: OPEN - unable to open file: ', Gname]);
    return
  end,


%  Put GRID NetCDF file into define mode.

  [status]=mexnc('redef',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['WRITE_MASK: REFDEF - unable to put into define mode.']);
    return
  end,

%  Define scope mask on RHO-points.

  if (define.Rscope),
    Var.name          = Vname.Rscope;
    Var.type          = ncdouble;
    Var.dimid         = [did.yr did.xr];
    Var.long_name     = 'adjoint sensitivity scope mask on RHO-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['inactive', blanks(1), ...
                         'active'];

    [varid,status]=nc_vdef(ncid,Var);
    clear Var
  end,

%  Define scope mask on U-points.

  if (define.Uscope),
    Var.name          = Vname.Uscope;
    Var.type          = ncdouble;
    Var.dimid         = [did.yu did.xu];
    Var.long_name     = 'adjoint sensitivity scope mask on U-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['inactive', blanks(1), ...
                         'active'];

    [varid,status]=nc_vdef(ncid,Var);
    clear Var
  end,

%  Define scope mask on V-points.

  if (define.Vscope),
    Var.name          = Vname.Vscope;
    Var.type          = ncdouble;
    Var.dimid         = [did.yv did.xv];
    Var.long_name     = 'adjoint sensitivity scope mask on V-points';
    Var.flag_values   = [0.0 1.0];
    Var.flag_meanings = ['inactive', blanks(1), ...
                         'active'];

    [varid,status]=nc_vdef(ncid,Var);
    clear Var
  end,

%  Leave definition mode.

  [status]=mexnc('enddef',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['WRITE_MASK: ENDDEF - unable to leave definition mode.']);
  end,

%  Close GRID NetCDF file.

  [status]=mexnc('close',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['WRITE_MASK: CLOSE - unable to close NetCDF file: ', Gname]);
  end,

end,

%---------------------------------------------------------------------------
%  Write out mask data into GRID NetCDF file.
%---------------------------------------------------------------------------

[status]=nc_write(Gname,Vname.Rscope,Rscope);
[status]=nc_write(Gname,Vname.Uscope,Uscope);
[status]=nc_write(Gname,Vname.Vscope,Vscope);

return
