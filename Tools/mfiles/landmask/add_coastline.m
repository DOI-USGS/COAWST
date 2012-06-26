function status=add_coastline(Gname, Clon, Clat);

%
% ADD_COASTLINE:  Adds coastline data to a ROMS Grid NetCDF file
%
% status=add_coastline(Gname, Clon, Clat)
%
% This add coasline data to GRID NetCDF file. This is the coastline
% data use for curvilinear grid and Land/Sea masking.
%
% On Input:
%
%    Gname       GRID NetCDF file name (character string).
%    Clon        Coastline longitude (degree_east).
%    Clat        Coastline latitude (degree_north).
%
% On Output:
%
%    status      Error flag.
%

% svn $Id: add_coastline.m 485 2010-07-07 18:10:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%


%---------------------------------------------------------------------------
% Inquire grid NetCDF file about coastline variables.
%---------------------------------------------------------------------------

got.Clon=0;  define.Clon=1;  Vname.Clon='lon_coast';
got.Clon=0;  define.Clat=1;  Vname.Clat='lat_coast';

[varnam,nvars]=nc_vname(Gname);
for n=1:nvars,
  name=deblank(varnam(n,:));
  switch name
    case {Vname.Clon}
      got.Clon=1;
      define.Clon=0;
    case {Vname.Clat}
      got.Clat=1;
      define.Clat=0;
  end,
end,

%---------------------------------------------------------------------------
%  If appropriate, define Land/Sea mask variables.
%---------------------------------------------------------------------------

if (define.Clon | define.Clat),

%  Inquire about dimensions.

  Dname.Clon='coast';
  Dname.Clat='coast';
  got.coast=0;

  [Dnames,Dsizes]=nc_dim(Gname);
  ndims=length(Dsizes);
  for n=1:ndims,
    dimid=n;
    name=deblank(Dnames(n,:));
    switch name
      case {Dname.Clon}
        Dsize.Clon=Dsizes(n);
        did.Clon=dimid;
        got.coast=1;
      case {Dname.Clat}
        Dsize.Clat=Dsizes(n);
        did.Clat=dimid;
        got.coast=1;
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
    error(['ADD_COASTLINE: OPEN - unable to open file: ', Gname]);
    return
  end,

%  Put GRID NetCDF file into define mode.

  [status]=mexnc('redef',ncid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['ADD_COASTLINE: REFDEF - unable to put into define mode.']);
    return
  end,

%  Define "coast" dimension.

   if (~got.coast),
     [did.Clon,status]=mexnc('def_dim',ncid,Dname.Clon,length(Clon));
     if (status ~= 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['ADD_COASTLINE: DEF_DIM - unable to define dimension: ',...
              Dname.Clon])
     end,
     did.Clat=did.Clon;
   end,

%  Define coastline longitude.

  if (define.Clon),
    Var.name      =  Vname.Clon;
    Var.type      =  ncdouble;
    Var.dimid     =  did.Clon;
    Var.long_name = 'Coastline longitude';
    Var.units     = 'degree_east';

    [varid,status]=nc_vdef(ncid,Var);
    clear Var
  end,

  if (define.Clat),
    Var.name      = Vname.Clat;
    Var.type      = ncdouble;
    Var.dimid     = did.Clat;
    Var.long_name = 'Coastline latitude';
    Var.units     = 'degree_north';

    [varid,status]=nc_vdef(ncid,Var);
    clear Var
  end,

%  Leave definition mode.

  [status]=mexnc('enddef',ncid);
  if (status ~= 0),
    error(['ADD_COASTLINE: ENDDEF - unable to leave definition mode.']);
  end,

%  Close GRID NetCDF file.

  [status]=mexnc('close',ncid);
  if (status ~= 0),
    error(['ADD_COASTLINE: CLOSE - unable to close NetCDF file: ', Gname]);
  end,

end,

%---------------------------------------------------------------------------
%  Write out coastline data into GRID NetCDF file.
%---------------------------------------------------------------------------

[status]=nc_write(Gname,Vname.Clon,Clon);
[status]=nc_write(Gname,Vname.Clat,Clat);

return
