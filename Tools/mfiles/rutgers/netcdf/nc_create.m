function ncid = nc_create(ncfile,mode,S)

%
% NC_CREATE:  It creates a New netCDF file
%
% ncid = nc_create(ncfile,mode,S)
%
% This function creates a new NetCDF file according to the file
% creation mode. If the structure S is provided, it defines the
% dimensions, global attributes, variables and attributes stored
% in the structure.  This structure can be created using:
%
%    S = nc_inq('roms_his.nc')
% or
%    S = ncinfo('roms_his.nc')        native Matlab function
%
% If the user wants different dimension values than those returned by
% the "nc_inq" or "ncinfo", see S.Dimensions(:).Length, make sure
% that such values are overwritten before calling this function.
% Therefore, it is possible to replicate the schema of a particular
% NetCDF file with identical or different dimensions. This facilitates
% the extraction of smaller subdomain data from a large (coarser)
% dataset.
%
% Warning: This function does no write the variable data. It just
%          creates the variables.
%
% On Input:
%
%    ncfile     NetCDF file name to create (string)
%
%    mode       NetCDF file creation mode (string):
%
%                 'clobber'       - overwrite existing files
%                 'noclobber'     - do not overwrite existing files
%                 'share'         - allow for synchronous file updates
%                 '64bit_offset'  - allow the creation of 64-bit files
%                                   instead of the classic format
%                 'NETCDF4'       - create a netCDF-4/HDF5 file
%                 'classic_model' - enforce classic model, has no effect
%                                   unless used in a bitwise-or with
%                                   'netcdf4'
%
%    S          NetCDF file Schema Structure (struct array)
%
% On Output:
%
%    ncid       NetCDF file ID. If no function output arguments are
%                 given, the file is closed after it is created.
%

% svn $Id: nc_create.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

if (nargin > 2),
  define = true;
else
  define = false;
end

%--------------------------------------------------------------------------
% Create NetCDF file.
%--------------------------------------------------------------------------

ncid = netcdf.create(ncfile,mode);

%--------------------------------------------------------------------------
% Define dimensions and variables available in dimension structure.
%--------------------------------------------------------------------------

if (define),

% Define dimensions.

  if (isfield(S,'Dimensions')),
    ndims = length(S.Dimensions);
  else
    disp(S);
    error([' NC_CREATE: unable to find ''Dimensions'' field',           ...
           ' in input structure, S.']);
  end

  for n=1:ndims,
    dname = char(S.Dimensions(n).Name);
    if (S.Dimensions(n).Unlimited),
      dlen = netcdf.getConstant('UNLIMITED');
    else
      dlen = S.Dimensions(n).Length;
    end
    Did.(dname) = netcdf.defDim(ncid,dname,dlen);
  end

% Define global attributes.

  if (isfield(S,'Attributes')),
    natts = length(S.Attributes);
  else
    disp(S);
    error([' NC_CREATE: unable to find ''Attributes'' field',           ...
           ' in input structure, S.']);
  end

  if (natts > 0)
    for n=1:natts,
      aname  = char(S.Attributes(n).Name);
      avalue = S.Attributes(n).Value;
      varid  = netcdf.getConstant('nc_global');
      netcdf.putAtt(ncid,varid,aname,avalue);  
    end
  end

% Define variables and their attributes.

  if (isfield(S,'Variables')),
    nvars = length(S.Variables);
  else
    disp(S);
    error([' NC_CREATE: unable to find ''Variables'' field',            ...
           ' in input structure, S.']);
  end
  
  got_nctype = isfield(S.Variables,'ncType');

  for n=1:nvars,
    nvdims = length(S.Variables(n).Dimensions);
    if (nvdims > 0)
      for i=1:nvdims,
        dname = char(S.Variables(n).Dimensions(i).Name);
        dimids(i) = Did.(dname);
      end
    else
      dimids=[];
    end
        
    if (~got_nctype),
      xtype = char(S.Variables(n).Datatype);
      switch (xtype)
        case 'int8'
          vtype = netcdf.getConstant('nc_byte');
        case 'uint8'
          vtype = netcdf.getConstant('nc_ubyte');
        case 'char'
          vtype = netcdf.getConstant('nc_char');
        case 'int16'
          vtype = netcdf.getConstant('nc_short');
        case 'uint16'
          vtype = netcdf.getConstant('nc_ushort');
        case 'int32'
          vtype = netcdf.getConstant('nc_int');
        case 'uint32'
          vtype = netcdf.getConstant('nc_uint');
        case 'single'
          vtype = netcdf.getConstant('nc_float');
        case 'double'
          vtype = netcdf.getConstant('nc_double');
        case 'int64'
          vtype = netcdf.getConstant('nc_int64');
        case 'uint64'
          vtype = netcdf.getConstant('nc_uint64');
        otherwise
          vtype = [];
      end
    else
      vtype = S.Variables(n).ncType;
    end  
  
    vname = char(S.Variables(n).Name);
    varid = netcdf.defVar(ncid,vname,vtype,dimids);

    nvatts = length(S.Variables(n).Attributes);
    if (nvatts > 0)
      for i=1:nvatts,
        aname  = char(S.Variables(n).Attributes(i).Name);
        avalue = S.Variables(n).Attributes(i).Value;
        netcdf.putAtt(ncid,varid,aname,avalue);  
      end
    end
    clear aname avalue dimids
  end

% End NetCDF define mode.

  netcdf.endDef(ncid);

end

%--------------------------------------------------------------------------
% Close NetCDF file.
%--------------------------------------------------------------------------

if (nargout < 1),
  netcdf.close(ncid);
end

return
