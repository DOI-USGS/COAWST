function V = nc_vinfo(ncfile, Vname)

%
% NC_VINFO:  Inquire information about requested NetCDF variable
%
% V = nc_vinfo(ncfile, Vname)
%
% This function gets information about requested NetCDF variable.
%
% On Input:
%
%    ncfile     NetCDF file name or URL name (character string)
%    Vname      Field variable name (character string)
%
% On Output:
%
%    V          Requested variable information (struct array):
%
%                 V.Filename    NetCDF file name (string)
%                 V.Name        variable names (string)
%                 V.Dimensions  variable dimensions (struct array):
%                                 V.Dimensions(:).Name
%                                 V.Dimensions(:).Length
%                                 V.Dimensions(:).Unlimited
%                 V.Size        variable size (double array)
%                 V.Attributes  variable attributes (struct array):
%                                 V.Attributes(:).Name
%                                 V.Attributes(:).Value
%                 V.Cgridtype   stagged C-grid type (struct array):
%                                 V.Cgridtype.Name
%                                 V.Cgridtype.Value
%                 V.Datatype    original variable data type (string)
%                 V.ncType      NetCDF data type
%

% svn $Id: nc_vinfo.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

if nargin < 2,
  error('NC_VINFO: variable name not provided ...');
end

% Choose NetCDF file interface.

[method,~,~] = nc_interface(ncfile);

switch(method),
  case {'native'}
    V = nc_vinfo_matlab(ncfile,Vname);           % Matlab native interface
  case {'java'}
    V = nc_vinfo_java  (ncfile,Vname);           % SNCTOOLS JAVA interface
  case {'mexnc'}
    V = nc_vinfo_mexnc (ncfile,Vname);           % MEXNC inteface
  otherwise
    error('NC_VINFO: unable to determine NetCDF processing interface');
end

return

%--------------------------------------------------------------------------

function V = nc_vinfo_matlab(ncfile, Vname)

%
% NC_VINFO_MATLAB:  Inquire information about requested NetCDF variable
%
% V = nc_vinfo_matlab(ncfile, Vname)
%
% This function gets information about requested variable from
% regular or OpenDAP URL NetCDF file. It uses native matlab function
% "ncinfo" that it is available in version 2012a or higher.
%

% Initialize output structure.

Attributes = struct('Name', '',                                         ...
                    'Value', []);
Cgridtype  = struct('Name', '',                                         ...
                    'Value', 0);
Dimensions = struct('Name', '',                                         ...
                    'Length', [],                                       ...
                    'Unlimited', false);

V = struct('Filename', '',                                              ...
           'Name', '',                                                  ...
           'Dimensions', Dimensions,                                    ...
           'Size', [],                                                  ...
           'Attributes', Attributes,                                    ...
           'Cgridtype', Cgridtype,                                      ...
           'Datatype', '',                                              ...
           'ncType', []);

% Inquire information about requested variable.

Info = ncinfo(ncfile,Vname);

% Determine C-grid type.

nvdims = length(Info.Dimensions);

Ctype.Name  = 'none';
Ctype.Value = 0;

wgrid = 0;

if (nvdims > 0),
  for n=1:nvdims,
    dname = Info.Dimensions(n).Name;
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

% Extract requested variable information.

V.Filename   = ncfile;
V.Name       = Vname;
V.Dimensions = Info.Dimensions;
V.Size       = Info.Size;
V.Attributes = Info.Attributes;
V.Cgridtype  = Ctype;
V.Datatype   = Info.Datatype;

% NetCDF data type.

vtype = char(Info.Datatype);
switch (vtype)
  case 'int8'
    V.ncType = netcdf.getConstant('nc_byte');
  case 'uint8'
    V.ncType = netcdf.getConstant('nc_ubyte');
  case 'char'
    V.ncType = netcdf.getConstant('nc_char');
  case 'int16'
    V.ncType = netcdf.getConstant('nc_short');
  case 'uint16'
    V.ncType = netcdf.getConstant('nc_ushort');
  case 'int32'
    V.ncType = netcdf.getConstant('nc_int');
  case 'uint32'
    V.ncType = netcdf.getConstant('nc_uint');
  case 'single'
    V.ncType = netcdf.getConstant('nc_float');
  case 'double'
    V.ncType = netcdf.getConstant('nc_double');
  case 'int64'
    V.ncType = netcdf.getConstant('nc_int64');
  case 'uint64'
    V.ncType = netcdf.getConstant('nc_uint64');
  otherwise
    V.ncType = [];
end

return

%--------------------------------------------------------------------------

function V = nc_vinfo_java(ncfile, Vname)

%
% NC_VINFO_JAVA:  Inquire information about requested NetCDF variable
%
% V = nc_vinfo_java(ncfile,Vname)
%
% This function gets information about requested variable from
% URL OpenDAP NetCDF file(s). It uses SNCTOOLS function "nc_info".
%

% Initialize output structure.

Attributes = struct('Name', '',                                         ...
                    'Value', []);
Cgridtype  = struct('Name', '',                                         ...
                    'Value', 0);
Dimensions = struct('Name', '',                                         ...
                    'Length', [],                                       ...
                    'Unlimited', false);

V = struct('Filename', '',                                              ...
           'Name', '',                                                  ...
           'Dimensions', Dimensions,                                    ...
           'Size', [],                                                  ...
           'Attributes', Attributes,                                    ...
           'Cgridtype', Cgridtype,                                      ...
           'Datatype', '',                                              ...
           'ncType', []);

%  Inquire information from URL NetCDF file, use SNCTOOLS.

Info = nc_varinfo(ncfile,Vname); 

% Determine C-grid type.

nvdims = length(Info.Dimension);

Ctype.Name  = 'none';
Ctype.Value = 0;

wgrid = 0;

if (nvdims > 0),
  for n=1:nvdims,
    dname = char(Info.Dimension{n});
    Unlimited(n) = nc_getdiminfo(ncfile,dname,'Unlimited');
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

% Extract requested variable information.

V.Filename   = ncfile;
V.Name       = Vname;
for n=1:nvdims
  V.Dimensions(n).Name      = Info.Dimension{n};
  V.Dimensions(n).Length    = Info.Size(n);
  V.Dimensions(n).Unlimited = Unlimited(n);
end
V.Size       = Info.Size;
V.Attributes = Info.Attribute;
V.Cgridtype  = Ctype;
V.Datatype   = Info.Datatype;

% NetCDF data type.

vtype = char(Info.Datatype);
switch (vtype)
  case 'int8'
    V.ncType = nc_constant('nc_byte');
  case 'uint8'
    V.ncType = nc_constant('nc_ubyte');
  case 'char'
    V.ncType = nc_constant('nc_char');
  case 'int16'
    V.ncType = nc_constant('nc_short');
  case 'uint16'
    V.ncType = nc_constant('nc_ushort');
  case 'int32'
    V.ncType = nc_constant('nc_int');
  case 'uint32'
    V.ncType = nc_constant('nc_uint');
  case 'single'
    V.ncType = nc_constant('nc_float');
  case 'double'
    V.ncType = nc_constant('nc_double');
  case 'int64'
    V.ncType = nc_constant('nc_int64');
  case 'uint64'
    V.ncType = nc_constant('nc_uint64');
  otherwise
    V.ncType = [];
end

return

%--------------------------------------------------------------------------

function V = nc_vinfo_mexnc(ncfile, Vname)

%
% NC_VINFO_MEXNC:  Inquire information about requested NetCDF variable
%
% V = nc_vinfo_mexnc(ncfile, Vname)
%
% This function gets information about requested variable from a
% NetCDF file. It uses MEXNC functions. Therefore, it cannot process
% an URL OpenDAP file.
%

% Initialize output structure.

Attributes = struct('Name', '',                                         ...
                    'Value', []);
Cgridtype  = struct('Name', '',                                         ...
                    'Value', 0);
Dimensions = struct('Name', '',                                         ...
                    'Length', [],                                       ...
                    'Unlimited', false);

V = struct('Filename', '',                                              ...
           'Name', '',                                                  ...
           'Dimensions', Dimensions,                                    ...
           'Size', [],                                                  ...
           'Attributes', Attributes,                                    ...
           'Cgridtype', Cgridtype,                                      ...
           'Datatype', '',                                              ...
           'ncType', []);

% Open NetCDF file.
 
[ncid]=mexnc('ncopen',ncfile,'nc_nowrite');
if (ncid == -1),
  error(['NC_VINFO: ncopen - unable to open file: ', ncfile]);
end
 
%  Supress all error messages from NetCDF.
 
mexnc('setopts',0);

%  Get variable ID.

[varid]=mexnc('ncvarid',ncid,Vname);
if (varid < 0),
  mexnc('ncclose',ncid);
  nc_inq(ncfile);
  disp('  ');
  error(['NC_VINFO: ncvarid - cannot find variable: ',Vname]);
end,

% Inquire about unlimmited dimension.

[~,~,~,recdim,status]=mexnc('ncinquire',ncid);
if (status == -1),
  error(['NC_VINFO: ncinquire - cannot inquire file: ',ncfile]);
end,
 
% Get information about requested variable.
 
[~,nctype,nvdims,dimids,nvatts,status]=mexnc('ncvarinq',ncid,varid);
if (status == -1),
  error(['NC_VINFO: ncvarinq - unable to inquire about variable: ',Vname]);
end,

% Inquire about dimensions.  Inverted the dimension to colunm major order.

V.Filename = ncfile;
V.Name     = Vname;

Ctype.Name  = 'none';
Ctype.Value = 0;

wgrid = 0;

if (nvdims > 0),
  n=0;
  for i=nvdims:-1:1
    n=n+1;
    [name,size,status]=mexnc('ncdiminq',ncid,dimids(i));
    if (status == -1),
      error(['NC_VINFO: ncdiminq - unable to inquire about dimension ', ...
             'ID: ',num2str(dimids(n))]);
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
      V.Dimensions(n).Name   = name;
      V.Dimensions(n).Length = size;
      if (dimids(n) == recdim),
         V.Dimensions(n).Unlimited = true;
      else
         V.Dimensions(n).Unlimited = false;
      end
      V.Size(n) = size;
    end
  end
else
  V.Dimensions = [];
  V.Size = 1;
end
if (wgrid == 5),
  Ctype.Name  = 'w-velocity point';
  Ctype.Value = 5;
end

% Inquire about the attributes.

if (nvatts > 0),
  for n=1:nvatts
    [attnam,status]=mexnc('inq_attname',ncid,varid,n-1);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_VINFO: inq_attname: error while inquiring attribute: ',...
             num2str(n)]);
    end 
    V.Attributes(n).Name = attnam;
    
    [atype,status]=mexnc('inq_atttype',ncid,varid,attnam);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_VINFO: inq_atttype - unable to inquire datatype for ', ...
             'attribute: ',attnam]);
    end

    switch (atype)
      case (nc_constant('nc_char'))
        [avalue,status]=mexnc('get_att_text',ncid,varid,attnam);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_VINFO: get_att_text - unable to read attribute ',  ...
                 '"',attnam,'"in variable: ',Vname,'.']);
        end
        V.Attributes(n).Value = avalue;  
      case (nc_constant('nc_int'))
        [avalue,status]=mexnc('get_att_int', ncid,varid,attnam);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_VINFO: get_att_int - unable to read attribute ',   ...
                 '"',attnam,'"in variable: ',Vname,'.']);
        end
        V.Attributes(n).Value = avalue;   
      case (nc_constant('nc_float'))
        [avalue,status]=mexnc('get_att_float',ncid,varid,attnam);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_VINFO: get_att_float - unable to read attribute ', ...
                 '"',attnam,'"in variable: ',Vname,'.']);
        end
        V.Attributes(n).Value = avalue;
      case (nc_constant('nc_double'))
        [avalue,status]=mexnc('get_att_double',ncid,varid,attnam);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_VINFO: get_att_double - unable to read attribute ',...
                 '"',attnam,'"in variable: ',Vname,'.']);
        end
        V.Attributes(n).Value = avalue;
    end
  end
else
  V.Attributes = [];
end

% Set data type.

V.Cgridtype  = Ctype;

switch(nctype)
  case nc_nat
    V.Datatype = '';
  case nc_byte
    V.Datatype = 'int8';
  case nc_ubyte
    V.Datatype = 'uint8';
  case nc_char
    V.Datatype = 'char';
  case nc_short
    V.Datatype = 'int16';
  case nc_ushort
    V.Datatype = 'uint16';
  case nc_int
    V.Datatype = 'int32';
  case nc_uint
    V.Datatype = 'uint32';
  case nc_float
    V.Datatype = 'single';
  case nc_double
    V.Datatype = 'double';
  case nc_int64
    V.Datatype = 'int64';
  case nc_uint64
    V.Datatype = 'uint64';
  otherwise
    V.Datatype = '';
end
V.ncType = nctype;

%   Close NetCDF file.

[status]=mexnc('ncclose',ncid);
if (status == -1),
  error('NC_VINFO: ncclose - unable to close NetCDF file.');
end

return

%--------------------------------------------------------------------------

function info = nc_varinfo(ncfile,varname,field)

%
% NC_VARINFO  Returns metadata about a specific NetCDF variable.
%
%   INFO = NC_VARINFO(NCFILE, VARNAME)
%
% Returns a metadata structure about the variable VARNAME in the netCDF
% file NCFILE.
%
%   INFO will have the following fields:
%
%     Name      - A string containing the name of the variable.
%     Datatype  - The datatype of the variable.
%     Unlimited - Either 1 if the variable has an unlimited dimension or 
%                   0 if not.
%     Dimension - a cell array with the names of the dimensions upon 
%                   which this variable depends.
%     Size      - Size of the variable.
%     Attribute - An array of structures corresponding to the attributes 
%                   defined for the specified variable.
%                         
%   INFO = NC_VARINFO(NCFILE,VARNAME,<'field'>)
%
% Returns only one of the above fields: Datatype, Unlimited, Dimension, 
% Size, Attribute. Handy for use in expressions.
%
% Each "Attribute" element is a struct itself and contains the following 
% fields:
%
%     Name      - A string containing the name of the attribute.
%     Datatype  - The datatype of the variable.
%     Value     - Value of the attribute.
%
% Adapted from SNCTOOLS:
%
%   https://mexcdf.svn.sourceforge.net/svnroot/mexcdf/snctools/trunk
%
%   This function was adapted from SNCTOOLS function "nc_getvarinfo" to
%   return the variable information with the original numerical native
%   precision.  The original function "nc_getvarinfo" unwisely converted
%   all the numerical values to double precision.  This is problematic
%   when dealing with the variable numerical attributes. Specially, when
%   using the NetCDF attributes:  "_FillValue", "missing_value",
%   "scale_factor", and "add_offset".  And to lesser extend when using
%   the attributes: "valid_min", "valid_max", and "valid_range".
%
%   This function it is self contained to avoid interactions with the
%   original SNCTOOLS functions. All function calls are private and
%   attached below.  Only the Java interface is used since it can also
%   process files on OpenDAP servers.  The original SNCTOOLS private
%   functions were renamed by removing the prefix "get".
%

info = nc_varinfo_java(ncfile,varname);

if nargin > 2,
  info = info.(field); 
end

return

function Dataset = nc_varinfo_java(ncfile,varname)

% Java backend for NC_VARINFO.

import ucar.nc2.dods.*     
import ucar.nc2.*         
                           
close_it = true;

% Try it as a local file.  If not a local file, try as via HTTP,
% then as dods.

if isa(ncfile,'ucar.nc2.NetcdfFile')
  jncid = ncfile;
  close_it = false;
elseif isa(ncfile,'ucar.nc2.dods.DODSNetcdfFile')
  jncid = ncfile;
  close_it = false;
elseif exist(ncfile,'file')
  fid = fopen(ncfile);
  ncfile = fopen(fid);
  fclose(fid);
  jncid = NetcdfFile.open(ncfile);
else
  try 
    jncid = NetcdfFile.open ( ncfile );
    catch %#ok<CTCH>
      try
        jncid = snc_opendap_open(ncfile);
        catch %#ok<CTCH>
          error (['NC_VARINFO: fileOpenFailure ', ...
                  'Could not open ''%s'' with java backend.' , ncfile]);
      end
  end
end

if isa(varname,'ucar.nc2.Variable')
  jvarid = varname;
else
  jvarid = jncid.findVariable(varname);
  if isempty(jvarid)
    close(jncid);
    error (['NC_VARINFO: badVariableName ', ...
            'Could not locate variable %s', varname]);
  end
end

% All the details are hidden here because we need the exact same
% functionality in "nc_info".

Dataset = nc_varidinfo_java(jvarid);

% If we were passed a java file id, don't close it upon exit.

if close_it
  close ( jncid );
end

return

function Dataset = nc_varidinfo_java(jvarid)

% NC_VARIDINFO_JAVA:  returns metadata structure for a netcdf variable
%
% This function is private to snctools.  It is called by nc_info and
% nc_varinfo, and uses the java API.
%
% USAGE:   Dataset = nc_varidinfo_java(jvarid);
% 
% PARAMETERS:
%
% Input:
%     jvarid:  
%         of type ucar.nc2.dods.DODSVariable
% Output:
%     Dataset:
%         struct of variable metadata

Attribute = struct('Name','','Nctype',0,'Datatype','','Value',NaN);

Dataset = struct('Name','','Nctype','','Datatype','','Unlimited',false, ...
                 'Dimension',{''},'Size',[],'Attribute',Attribute,      ...
                 'Chunking',[],'Shuffle',0,'Deflate',0);

Dataset.Name = char ( jvarid.getName() );

% Get the datatype, store as an integer.

datatype = char(jvarid.getDataType().toString());

switch ( datatype )
  case 'double'
    Dataset.Nctype = nc_constant('nc_double');
    Dataset.Datatype = datatype;
  case 'float'
    Dataset.Nctype = nc_constant('nc_float');
    Dataset.Datatype = 'single';
  case {'int','long'}
    Dataset.Nctype = nc_constant('nc_int');
    Dataset.Datatype = 'int32';
  case 'short'
    Dataset.Nctype = nc_constant('nc_short');
    Dataset.Datatype = 'int16';
  case 'String'
    Dataset.Nctype = 12;           % Apparently, DODSNetcdfFile returns
    Dataset.Datatype = 'string';   % 'String', while NetcdfFile returns 'char'
  case 'char'
    Dataset.Nctype = nc_constant('nc_char');
    Dataset.Datatype = 'char';
  case 'byte'
    Dataset.Nctype = nc_constant('nc_byte');
    Dataset.Datatype = 'int8';
  otherwise
    error (['NC_VARINFO: unhandledDatatype ', ...
            '%s:  unhandled datatype ''%s''\n', datatype]);
end

% Determine if it is unlimited or not.

Dataset.Unlimited = double ( jvarid.isUnlimited() );

% Retrieve the dimensions.

dims = jvarid.getDimensions();
nvdims = dims.size();
Dimension = cell(1,nvdims);

for j = 1:nvdims
  theDim = jvarid.getDimension(j-1);
  Dimension{j} = char ( theDim.getName() );
end
Dataset.Dimension = Dimension;

% Get the size of the variable

if nvdims == 0
  Dataset.Size = 1;
else
  Size = jvarid.getShape();
  Dataset.Size = Size';
end

if nc_getpref('PRESERVE_FVD')
  Dataset.Dimension = fliplr(Dataset.Dimension);
  Dataset.Size = fliplr(Dataset.Size);
end

% Get the list of attributes.

j_att_list = jvarid.getAttributes();
Dataset.Attribute = nc_attsinfo_java(j_att_list);

return

function Attribute = nc_attsinfo_java(j_att_list)

% NC_ATTSINFO_JAVA:  returns metadata about netcdf attributes
%
% USAGE:  Attribute = nc_attsinfo_java(j_att_list);
%
% PARAMETERS:
%
% Input:
%     j_att_list:
%         Of type "java.util.ArrayList".  Each list member is of type
%         "ucar.nc2.Attribute"
% Output:
%     Attribute:
%         Structure array of attribute metadata.  The fields are 
%         
%         Name
%         Nctype (backwards compatibility)
%         Datatype
%         Value

j_att_iterator = j_att_list.listIterator();
j = 0;

Attribute = struct('Name','','Nctype',0,'Datatype','','Value',[]);

while 1
    
% This throws an exception when we've reached the end of the list.

  try
    jatt = j_att_iterator.next();
    catch %#ok<CTCH>
      break;
  end
    
  j = j + 1;
  Attribute(j) = nc_attinfo_java(jatt);

end

if j == 0
  Attribute = [];
end

return

function Attribute = nc_attinfo_java(jatt)

% NC_ATTINFO_JAVA:  return metadata about netcdf attribute

Attribute = struct('Name','','Nctype',0,'Datatype','','Value',[]);
    
Attribute.Name = char(jatt.getName());
    
datatype = char(jatt.getDataType().toString());

switch ( datatype )
  case 'double'
    Attribute.Nctype = 6; %#ok<*AGROW>
    Attribute.Datatype = 'double';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = double(values)';
        
  case 'float'
    Attribute.Nctype = 5;
    Attribute.Datatype = 'single';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = single(values)';
        
  case 'String'
    Attribute.Nctype = 12;
    Attribute.Datatype = 'string';
    shape = double(jatt.getLength);
    Attribute.Value = snc_pp_strings( jatt, jatt.getValues(), shape) ;
        
  case 'char'
    Attribute.Nctype = 2;
    Attribute.Datatype = 'char';
    Attribute.Value = char ( jatt.getStringValue());
        
  case 'byte'
    Attribute.Nctype = 1;
    Attribute.Datatype = 'int8';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = int8(values)';
        
  case 'short'
    Attribute.Nctype = 3;
    Attribute.Datatype = 'int16';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = int16(values)';
        
  case 'int'
    Attribute.Nctype = 4;
    Attribute.Datatype = 'int32';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = int32(values)';
        
  case 'long'
    Attribute.Nctype = 4;
    Attribute.Datatype = 'int64';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = int64(values)';
        
  otherwise
    error(['NC_VARINFO: unhandledDatatype  ', ...
           'Unhandled attribute datatype ''%s''\n', datatype]);
end

return

function value = nc_getpref(name)

%   PREF = NC_GETPREF(NAME) returns the value of an SNCTOOLS preference.
%   
%   This routine should not be called directly.

persistent PRESERVE_FVD

if isempty(PRESERVE_FVD)
    PRESERVE_FVD = getpref('SNCTOOLS','PRESERVE_FVD',false);
end

if strcmp(name,'PRESERVE_FVD')
    value = PRESERVE_FVD;
else
    error('unrecognized input to NC_GETPREF');
end

function jncid = snc_opendap_open(ncfile)

% SNC_OPENDAP_OPEN Open a connection to an OPeNDAP URL.  If the URL is
% password-protected, we will have to coerce netcdf-java to supply the
% credentials.

import ucar.nc2.dods.*  

% Is it a username/password protected URL?

pat = '(?<protocol>https{0,1})://(?<username>[^:]+):(?<password>[^@]+)@(?<host>[^/]+)(?<relpath>.*)';
parts = regexp(ncfile,pat,'names');

if numel(parts) == 0
  jncid = DODSNetcdfFile(ncfile);
else           % SncCreds is a custom java class supplied with SNCTOOLS.
  credentials = SncCreds(parts.username,parts.password);
  client = ucar.nc2.util.net.HttpClientManager.init(credentials,'snctools');
  opendap.dap.DConnect2.setHttpClient(client);
  ucar.unidata.io.http.HTTPRandomAccessFile.setHttpClient(client);
  ucar.nc2.dataset.NetcdfDataset.setHttpClient(client);
    
  jncid = DODSNetcdfFile(ncfile);
end

return

function values = snc_pp_strings(jobj,jdata,shape)

% Post process NC_STRING data into cell arrays.

if isempty(jdata)
  values = {''};
  return;
elseif strcmp(version('-release'),'14')
  % In R14, we must use the 2.2.x release of java.  No access to the
  % "getObject" method.  Assuming a single-valued string.
  values = {char(jobj.getStringValue())};
  return;
end

% Java says that the variable is laid out in row-major order.

if numel(shape) == 1
  values = cell([1 shape]);
else
  values = cell(shape);
end

for j = 1:prod(shape)
  values{j} = jdata.getObject(j-1);
end

return

