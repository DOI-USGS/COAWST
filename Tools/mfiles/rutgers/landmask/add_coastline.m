function add_coastline(ncfile, Clon, Clat)

%
% ADD_COASTLINE:  Adds coastline data to a ROMS Grid NetCDF file
%
% add_coastline(ncfile, Clon, Clat)
%
% Adds coastline data to a ROMS Grid NetCDF file.  The coastline
% data is used to facilitate Land/Sea masking.
%
% On Input:
%
%    ncfile      GRID NetCDF file name (string)
%    Clon        Coastline longitude (1D array; degree_east)
%    Clat        Coastline latitude  (1D array; degree_north)
%
%

% svn $Id: add_coastline.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire grid NetCDF file about coastline variables.
%--------------------------------------------------------------------------

got.Clon=false;  define.Clon=true;  Vname.Clon='lon_coast';
got.Clon=false;  define.Clat=true;  Vname.Clat='lat_coast';

V=nc_vnames(ncfile);
nvars=length(V.Variables);
for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch name
    case {Vname.Clon}
      got.Clon=true;
      define.Clon=false;
    case {Vname.Clat}
      got.Clat=true;
      define.Clat=false;
  end
end

%--------------------------------------------------------------------------
%  If appropriate, define Land/Sea mask variables.
%--------------------------------------------------------------------------

if (define.Clon || define.Clat),

%  Inquire about dimensions.

  Dname.Clon='coast';
  Dname.Clat='coast';
  got.coast=false;

  D=nc_dinfo(ncfile);
  ndims=length(D);
  for n=1:ndims,
    name=char(D(n).Name);
    switch name
      case {Dname.Clon}
        Dsize.Clon=D(n).Name;
        did.Clon=D(n).dimid;
        got.coast=true;
      case {Dname.Clat}
        Dsize.Clat=D(n).Name;
        did.Clat=D(n).dimid;
        got.coast=true;
    end
  end

%  Open GRID NetCDF file and put into define mode.

  ncid=netcdf.open(ncfile,'WRITE');
  netcdf.reDef(ncid);

%  Define "coast" dimension.

  if (~got.coast),
    did.Clon=netcdf.defDim(ncid,Dname.Clon,length(Clon));
    did.Clat=did.Clon;
  end

%  Define coastline longitude.

  if (define.Clon),
    varid=netcdf.defVar(ncid,Vname.Clon,                                ...
                        netcdf.getConstant('nc_double'),did.Clon);
    netcdf.putAtt(ncid,varid,'long_name','coastline longitude');
    netcdf.putAtt(ncid,varid,'units','degree_east');
  end

  if (define.Clat),
    varid=netcdf.defVar(ncid,Vname.Clat,                                ...
                        netcdf.getConstant('nc_double'),did.Clat);
    netcdf.putAtt(ncid,varid,'long_name','coastline latitude');
    netcdf.putAtt(ncid,varid,'units','degree_north');
  end

%  Leave definition mode and close NetCDF file.

  netcdf.endDef(ncid);
  netcdf.close(ncid);

end

%--------------------------------------------------------------------------
%  Write out coastline data into GRID NetCDF file.
%--------------------------------------------------------------------------

ncwrite(ncfile,Vname.Clon,double(Clon));
ncwrite(ncfile,Vname.Clat,double(Clat));

return
