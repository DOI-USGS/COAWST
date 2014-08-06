function nc_append(ncfile,S)

%
% NC_APPEND:  It appends new variable(s) to a existing NetCDF file
%
% ncid = nc_append(ncfile,mode,S)
%
% This function defines/appends new variables to an existing NetCDF
% file.  The variable(s) metadata is provided in input structure
% S, which includes information about variable type, dimensions and
% attributes.  It checks if the required variable dimensions are
% available.  If not, it defines needed dimensions.
%
% The input structure needs to have the same global dimensions
% and variable fields (schema) as those returned by:
%
%    S = nc_inq(ncfile)
% or
%    S = ncinfo(ncfile)          native Matlab function
%
% Warning:  this function will define all the variables in input
%           structure S that are not present in existing NetCDF
% file.  You may build the structure by hand or subsample the
% structure returned by 'nc_inq' or 'ncinfo' to contain just the
% desired variables to define/append.
%
% This function does no write the variable data. It just define the
% variables.
%
% On Input:
%
%    ncfile     Existing NetCDF file name (string)
%
%    S          NetCDF file Schema Structure (struct array)
%
%

% svn $Id: nc_append.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

if (~isstruct(S)),
  error('NC_APPEND: input argument ''S'' is not a structure');
end

if (~isfield(S,'Dimensions')),
  disp(S);
  error(['NC_APPEND: unable to find ''Dimensions'' field',              ...
         ' in input structure, S.']);
end

if (~isfield(S,'Variables')),
  disp(S);
  error(['NC_APPEND: unable to find ''Variables'' field',              ...
         ' in input structure, S.']);
end

%--------------------------------------------------------------------------
% Get existing file information.
%--------------------------------------------------------------------------

I = nc_inq(ncfile);

%--------------------------------------------------------------------------
% Open input NetCDF and put it in define mode.
%--------------------------------------------------------------------------

ncid = netcdf.open(ncfile,'WRITE');

netcdf.reDef(ncid);

%--------------------------------------------------------------------------
% Define variable(s) not available in existing NetCDF file.
%--------------------------------------------------------------------------

nvars = length(S.Variables);

for n=1:nvars,

  vname = char(S.Variables(n).Name);

  if (~strcmp({I.Variables.Name}, vname)),
    nvdims = length(S.Variables(n).Dimensions);
    
% If applicable, define needed variable dimensions.

    for i=1:nvdims,
      dname = char(S.Variables(n).Dimensions(i).Name);
      if (~strcmp({I.Dimensions.Name}, dname)),
        if (S.Variables(n).Dimensions(i).Unlimited),
          dlen = netcdf.getConstant('UNLIMITED');
        else
          dlen = S.Variables(n).Dimensions(i).Length;
        end
        Did.(dname) = netcdf.defDim(ncid,dname,dlen);
      else
        Did.(dname) = netcdf.inqDimID(ncid,dname);
      end
    end

% Define variable.

    got_nctype = isfield(S.Variables,'ncType');

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
  
    varid = netcdf.defVar(ncid,vname,vtype,dimids);

    clear dimids
    
% Define variable attributes.
  
    nvatts = length(S.Variables(n).Attributes);
    if (nvatts > 0)
      for i=1:nvatts,
        aname  = char(S.Variables(n).Attributes(i).Name);
        avalue = S.Variables(n).Attributes(i).Value;
        netcdf.putAtt(ncid,varid,aname,avalue);  
      end
      clear aname avalue dimids
    end
  
  end
end

%--------------------------------------------------------------------------
% End NetCDF define mode and close NetCDF file.
%--------------------------------------------------------------------------

netcdf.endDef(ncid);
netcdf.close(ncid);

return
