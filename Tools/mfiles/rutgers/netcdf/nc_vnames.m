function S = nc_vnames(fname)

%
% NC_VNAMES:  Gets the names of all variables in a NetCDF file
%
% S = nc_vnames(fname)
%
% This function gets the names and informaton of all variables available
% in a NetCDF file.
%
% On Input:
%
%    fname      NetCDF file name or URL name (character string)
%
% On Output:
%
%    S          NetCDF file variables information (struct array):
%
%                 S.Filename    NetCDF file name (string)
%                 S.Variables   variables information (struct array):
%
%                                 S.Variables(:).Name
%
%                                 S.Variables(:).Dimensions(:).Name
%                                 S.Variables(:).Dimensions(:).Length
%                                 S.Variables(:).Dimensions(:).Unlimited
%
%                                 S.Variables(:).Size
%
%                                 S.Variables(:).Attributes(:).Name
%                                 S.Variables(:).Attributes(:).Value
%
%                                 S.Variables(:).Cgridtype.Name
%                                 S.Variables(:).Cgridtype.Value
%
%                                 S.Variables(:).Datatype
%
%                Notice that the number of variables in the structure is:
%
%                  nvars = length(S.Variables);
%

% svn $Id: nc_vnames.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Choose NetCDF file interface.

[method,~,~] = nc_interface(fname);

switch(method),
  case {'native'}
    S = nc_vnames_matlab(fname);                % Matlab native interface
  case {'java'}
    S = nc_vnames_java  (fname);                % SNCTOOLS JAVA interface
  case {'mexnc'}
    S = nc_vnames_mexnc (fname);                % MEXNC inteface
  otherwise
    error('NC_VNAMES: unable to determine NetCDF processing interface');
end

return

%--------------------------------------------------------------------------

function S = nc_vnames_matlab(fname)

%
% NC_VNAMES_MATLAB:  Get the names of all variables in a NetCDF file
%
% S = nc_vnames_matlab(fname)
%
% This function gets information about all the variable available in
% a regular or OpenDAP URL NetCDF file. It uses the native matlab
% function "ncinfo" that it is available in version 2012a or higher.
%

% Initialize output structure.

Attributes = struct('Name', '',                                         ...
                    'Value', []);
Cgridtype  = struct('Name', '',                                         ...
                    'Value', 0);
Dimensions = struct('Name', '',                                         ...
                    'Length', [],                                       ...
                    'Unlimited', false);
Variables  = struct('Name', '',                                         ...
                    'Dimensions', Dimensions,                           ...
                    'Size', [],                                         ...
                    'Attributes', Attributes,                           ...
                    'Cgridtype', Cgridtype,                             ...
                    'Datatype', '',                                     ...
                    'ncType', []);

S = struct('Filename', '',                                              ...
           'Variables', Variables);

% Inquire information about requested variable.

Info = ncinfo(fname); 

% Extract variables information.

nvars = length(Info.Variables);

S.Filename = fname;

for n=1:nvars,
  S.Variables(n).Name       = Info.Variables(n).Name;
  S.Variables(n).Dimensions = Info.Variables(n).Dimensions;
  S.Variables(n).Size       = Info.Variables(n).Size;
  S.Variables(n).Attributes = Info.Variables(n).Attributes;

% Determine C-grid type.

  nvdims = length(Info.Variables(n).Dimensions);

  Ctype.Name  = 'none';
  Ctype.Value = 0;

  wgrid = 0;

  if (nvdims > 0),
    for i=1:nvdims,
      dname = Info.Variables(n).Dimensions(i).Name;
      switch (dname)
        case {'xi_rho', 'eta_rho', 'lon_rho', 'lat_rho', 'lon', 'lat'}
          Ctype.Name  = 'density point';
          Ctype.Value = 1;
        case {'xi_psi', 'eta_psi', 'lon_psi', 'lat_psi'}
          Ctype.Name  = 'streamfunction point';
          Ctype.Value = 2;
        case {'xi_u', 'eta_u', 'lon_u', 'lat_u'}
          Ctype.Name  = 'u-velocity point';
          Ctype.Value = 3;
        case {'xi_v', 'eta_v', 'lon_v', 'lat_v'}
          Ctype.Name  = 'v-velocity point';
          Ctype.Value = 4;
        case 's_w'
          wgrid = 5;
      end
    end  
    if (wgrid == 5),
      Ctype.Name  = 'w-velocity point';
      Ctype.Value = 5;
    end
  end
  S.Variables(n).Cgridtype = Ctype;

  S.Variables(n).Datatype  = Info.Variables(n).Datatype;

% NetCDF data type.

  vtype = char(Info.Variables(n).Datatype);
  switch (vtype)
    case 'int8'
      S.Variables(n).ncType = netcdf.getConstant('nc_byte');
    case 'uint8'
      S.Variables(n).ncType = netcdf.getConstant('nc_ubyte');
    case 'char'
      S.Variables(n).ncType = netcdf.getConstant('nc_char');
    case 'int16'
      S.Variables(n).ncType = netcdf.getConstant('nc_short');
    case 'uint16'
      S.Variables(n).ncType = netcdf.getConstant('nc_ushort');
    case 'int32'
      S.Variables(n).ncType = netcdf.getConstant('nc_int');
    case 'uint32'
      S.Variables(n).ncType = netcdf.getConstant('nc_uint');
    case 'single'
      S.Variables(n).ncType = netcdf.getConstant('nc_float');
    case 'double'
      S.Variables(n).ncType = netcdf.getConstant('nc_double');
    case 'int64'
      S.Variables(n).ncType = netcdf.getConstant('nc_int64');
    case 'uint64'
      S.Variables(n).ncType = netcdf.getConstant('nc_uint64');
    otherwise
      S.Variables(n).ncType = [];
  end
end

return

%--------------------------------------------------------------------------

function S = nc_vnames_java(fname)

%
% NC_VNAMES_JAVA:  Get the names of all variables in a NetCDF file
%
% S = nc_vname_java(fname)
%
% This function gets information about all the variable available in
% an OpenDAP URL NetCDF file(s).
%

% Initialize output structure.

Attributes = struct('Name', '',                                         ...
                    'Value', []);
Cgridtype  = struct('Name', '',                                         ...
                    'Value', 0);
Dimensions = struct('Name', '',                                         ...
                    'Length', [],                                       ...
                    'Unlimited', false);
Variables  = struct('Name', '',                                         ...
                    'Dimensions', Dimensions,                           ...
                    'Size', [],                                         ...
                    'Attributes', Attributes,                           ...
                    'Cgridtype', Cgridtype,                             ...
                    'Datatype', '',                                     ...
                    'ncType', []);

S = struct('Filename', '',                                              ...
           'Variables', Variables);

% Inquire information from URL NetCDF file, use SNCTOOLS.

Info = nc_info(fname); 

% Extract requested variable information.

nvars = length(Info.Dataset);

S.Filename = fname;

for n=1:nvars,
  S.Variables(n).Name = Info.Dataset(n).Name;

  nvdims = length(Info.Dataset(n).Dimension);
  for i=1:nvdims,
    S.Variables(n).Dimensions(i).Name   = Info.Dataset(n).Dimension;
    S.Variables(n).Dimensions(i).Length = Info.Dataset(n).Size;
  end

  S.Variables(n).Size = Info.Dataset(n).Size;

  nvatts = length(Info.Dataset(n).Attribute);
  for i=1:nvatts,
    S.Variables(n).Attributes(i).Name =Info.Dataset(n).Attribute(i).Name;
    S.Variables(n).Attributes(i).Value=Info.Dataset(n).Attribute(i).Value;
  end
  
% Determine C-grid type.

  Ctype.Name  = 'none';
  Ctype.Value = 0;
  
  wgrid = 0;

  if (nvdims > 0),
    for i=1:nvdims,
      dname = char(Info.Dataset(n).Dimension{i});
      Unlimited = nc_getdiminfo(fname,dname,'Unlimited');
      S.Variables(n).Dimensions(i).Unlimited = Unlimited;
      switch (dname)
        case {'xi_rho', 'eta_rho', 'lon_rho', 'lat_rho', 'lon', 'lat'}
          Ctype.Name  = 'density point';
          Ctype.Value = 1;
        case {'xi_psi', 'eta_psi', 'lon_psi', 'lat_psi'}
          Ctype.Name  = 'streamfunction point';
          Ctype.Value = 2;
        case {'xi_u', 'eta_u', 'lon_u', 'lat_u'}
          Ctype.Name  = 'u-velocity point';
          Ctype.Value = 3;
        case {'xi_v', 'eta_v', 'lon_v', 'lat_v'}
          Ctype.Name  = 'v-velocity point';
          Ctype.Value = 4;
        case 's_w'
          wgrid = 5;
      end
    end  
    if (wgrid == 5),
      Ctype.Name  = 'w-velocity point';
      Ctype.Value = 5;
    end
  end
  S.Variables(n).Cgridtype = Ctype;

  S.Variables(n).Datatype  = Info.Dataset(n).Datatype;

% NetCDF data type.

  vtype = char(Info.Dataset(n).Datatype);
  switch (vtype)
    case 'int8'
      S.Variables(n).ncType = nc_constant('nc_byte');
    case 'uint8'
      S.Variables(n).ncType = nc_constant('nc_ubyte');
    case 'char'
      S.Variables(n).ncType = nc_constant('nc_char');
    case 'int16'
      S.Variables(n).ncType = nc_constant('nc_short');
    case 'uint16'
      S.Variables(n).ncType = nc_constant('nc_ushort');
    case 'int32'
      S.Variables(n).ncType = nc_constant('nc_int');
    case 'uint32'
      S.Variables(n).ncType = nc_constant('nc_uint');
    case 'single'
      S.Variables(n).ncType = nc_constant('nc_float');
    case 'double'
      S.Variables(n).ncType = nc_constant('nc_double');
    case 'int64'
      S.Variables(n).ncType = nc_constant('nc_int64');
    case 'uint64'
      S.Variables(n).ncType = nc_constant('nc_uint64');
    otherwise
      S.Variables(n).ncType = [];
  end
end

return

%--------------------------------------------------------------------------

function S = nc_vnames_mexnc(fname)

%
% NC_VNAMES_MEXCDF:  Get the names of all variables in a NetCDF file
%
% S = nc_vname_mexcdf(fname)
%
% This function gets information about all the variable available in
% a NetCDF file. It cannot process URL OpenDAP file(s).

% Initialize output structure.

Attributes = struct('Name', '',                                         ...
                    'Value', []);
Cgridtype  = struct('Name', '',                                         ...
                    'Value', 0);
Dimensions = struct('Name', '',                                         ...
                    'Length', [],                                       ...
                    'Unlimited', false);
Variables  = struct('Name', '',                                         ...
                    'Dimensions', Dimensions,                           ...
                    'Size', [],                                         ...
                    'Attributes', Attributes,                           ...
                    'Cgridtype', Cgridtype,                             ...
                    'Datatype', '',                                     ...
                    'ncType', []);

S = struct('Filename', '',                                              ...
           'Variables', Variables);

% Open NetCDF file.
 
[ncid]=mexnc('ncopen',fname,'nc_nowrite');
if (ncid == -1),
  error(['NC_VNAMES: ncopen - unable to open file: ', fname]);
end
 
% Supress all error messages from NetCDF.
 
mexnc('setopts',0);

% Inquire about the contents of NetCDf file. Display information.

[~,nvars,~,recdim,status]=mexnc('ncinquire',ncid);
if (status == -1),
  error([ 'NC_VNAMES: ncinquire - error while inquiring file: ' fname]);
end,

% Get information about all variables.

S.Filename = fname;

for n=1:nvars
  varid=n-1;
  [name,nctype,nvdims,dimids,nvatts,status]=mexnc('ncvarinq',ncid,varid);
  if (status == -1),
    error(['NC_VNAMES: ncvarinq - unable to inquire about variable: ', ...
           vname]);
  end,
  S.Variables(n).Name = name;
  
% Inquire about dimensions.  Inverted the dimension to colunm major order.

  Ctype.Name  = 'none';
  Ctype.Value = 0;

  wgrid = 0;

  if (nvdims > 0),
    m=0;
    for i=nvdims:-1:1
      m=m+1;
      [name,size,status]=mexnc('ncdiminq',ncid,dimids(i));
      if (status == -1),
        error(['NC_VNAMES: ncdiminq - unable to inquire about ',        ...
               'dimension ID: ',num2str(dimids(m))]);
      else
        switch (name)
          case {'xi_rho', 'eta_rho', 'lon_rho', 'lat_rho', 'lon', 'lat'}
            Ctype.Name  = 'density point';
            Ctype.Value = 1;
          case {'xi_psi', 'eta_psi', 'lon_psi', 'lat_psi'}
            Ctype.Name  = 'streamfunction point';
            Ctype.Value = 2;
          case {'xi_u', 'eta_u', 'lon_u', 'lat_u'}
            Ctype.Name  = 'u-velocity point';
            Ctype.Value = 3;
          case {'xi_v', 'eta_v', 'lon_v', 'lat_v'}
            Ctype.Name  = 'v-velocity point';
            Ctype.Value = 4;	 
          case 's_w'
            wgrid = 5;
	end
        S.Variables(n).Dimensions(m).Name   = name;
        S.Variables(n).Dimensions(m).Length = size;
        if (dimids(m) == recdim),
          S.Variables(n).Dimensions(m).Unlimited = true;
        else
          S.Variables(n).Dimensions(m).Unlimited = false;
        end
        S.Variables(n).Size(m) = size;
      end
    end
  else
    S.Variables(n).Dimensions = [];
    S.Variables(n).Size = 1;
  end
  if (wgrid == 5),
    Ctype.Name  = 'w-velocity point';
    Ctype.Value = 5;
  end

% Inquire about the attributes.

  if (nvatts > 0),
    for m=1:nvatts
      [attnam,status]=mexnc('inq_attname',ncid,varid,m-1);
      if (status < 0),
        disp('  ');
        disp(mexnc('strerror',status));
        error(['NC_VNAMES: inq_attname: error while inquiring ',        ...
               'attribute: ',num2str(m)]);
      end 
      S.Variables(n).Attributes(m).Name = attnam;
    
      [atype,status]=mexnc('inq_atttype',ncid,varid,attnam);
      if (status < 0),
        disp('  ');
        disp(mexnc('strerror',status));
        error(['NC_VNAMES: inq_atttype - unable to inquire datatype ',  ...
               'for attribute: ',attnam]);
      end

      switch (atype)
        case (nc_constant('nc_char'))
          [avalue,status]=mexnc('get_att_text',ncid,varid,attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_VNAMES: get_att_text - unable to read ',         ...
                   'attribute "',attnam,'"in variable: ',vname,'.']);
          end
          S.Variables(n).Attributes(m).Value = avalue;  
        case (nc_constant('nc_int'))
          [avalue,status]=mexnc('get_att_int', ncid,varid,attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_VNAMES: get_att_int - unable to read ',          ...
                   'attribute "',attnam,'"in variable: ',vname,'.']);
          end
          S.Variables(n).Attributes(m).Value = avalue;   
        case (nc_constant('nc_float'))
          [avalue,status]=mexnc('get_att_float',ncid,varid,attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_VINFO: get_att_float - unable to read ',         ...
                   'attribute "',attnam,'"in variable: ',vname,'.']);
          end
          S.Variables(n).Attributes(m).Value = avalue;
        case (nc_constant('nc_double'))
          [avalue,status]=mexnc('get_att_double',ncid,varid,attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_VINFO: get_att_double - unable to read ',        ...
                   'attribute "',attnam,'"in variable: ',vname,'.']);
          end
          S.Variables(n).Attributes(m).Value = avalue;
      end
    end
  else
    S.Variables(n).Attributes = [];
  end

% Set data type.

  S.Variables(n).Cgridtype  = Ctype;

  switch(nctype)
    case nc_nat
      S.Variables(n).Datatype = '';
    case nc_byte
      S.Variables(n).Datatype = 'int8';
    case nc_ubyte
      S.Variables(n).Datatype = 'uint8';
    case nc_char
      S.Variables(n).Datatype = 'char';
    case nc_short
      S.Variables(n).Datatype = 'int16';
    case nc_ushort
      S.Variables(n).Datatype = 'uint16';
    case nc_int
      S.Variables(n).Datatype = 'int32';
    case nc_uint
      S.Variables(n).Datatype = 'uint32';
    case nc_float
      S.Variables(n).Datatype = 'single';
    case nc_double
      S.Variables(n).Datatype = 'double';
    case nc_int64
      S.Variables(n).Datatype = 'int64';
    case nc_uint64
      S.Variables(n).Datatype = 'uint64';
    otherwise
      S.Variables(n).Datatype = '';
  end

  S.Variables(n).ncType = nctype;
end

%   Close NetCDF file.

[status]=mexnc('ncclose',ncid);
if (status == -1),
  error('NC_VNAMES: ncclose - unable to close NetCDF file.');
end

return
