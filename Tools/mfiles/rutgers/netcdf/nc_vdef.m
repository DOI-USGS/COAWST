function [varid, status] = nc_vdef(ncid, Var)

%
% NC_VDEF:   Create a ROMS variable in a NetCDF file
%
% [varid, status] = nc_vdef(ncid, Var)
%
% This function defines a variable into NetCDF file.
%
% On Input:
%
%    ncid        NetCDF file ID
%
%    Var         Variable information (structure array):
%
%                  Var.name  => variable name (string)
%
%                  Var.type  => external data type (string or numeric):
%                               'byte'   or  nc_byte
%                               'char'   or  nc_char
%                               'short'  or  nc_short
%                               'int'    or  nc_int
%                               'float'  or  nc_float
%                               'double' or  nc_double
%                                     
%                  Var.dimid => dimension IDs (numeric), if not present
%                               or empthy [] the variable is a scalar
%                  
%                  Any other field in the structure, if any,
%                  is proccessed as a variable attributes. For
%                  example:
%
%                  Var.long_name  => "long_name" attribute (string)
%                  Var.FillValue  => "_FillValue" attribute (number)
%
%                  the values of the attribute can be numeric (scalar
%                  or vector) or characters (array or cell array).
%
% On Output:
%
%    varid       Variable ID
%
%    status      Error flag
%

% svn $Id: nc_vdef.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Choose NetCDF file interface.

[method,~,~] = nc_interface();

switch(method),
  case {'native'}
    [varid, status] = nc_vdef_matlab(ncid, Var);
  case {'mexnc'}
    [varid, status] = nc_vdef_mexnc (ncid, Var);
  case {'java'}
    error('NC_VDEF: it is not possible to write into an OpenDAP file.');
  otherwise
    error('NC_VDEF: unable to determine NetCDF processing interface.');
end

return

%--------------------------------------------------------------------------

function [varid, status] = nc_vdef_matlab(ncid, Var)

%
% NC_VDEF_MEXNC:   Create a ROMS variable in a NetCDF file
%
% [varid, status] = nc_vdef(ncid, Var)
%
% This function defines a variable into NetCDF file using the native
% Matlab inteface.
%
% On Input:
%
%    ncid        NetCDF file ID
%
%    Var         Variable information (structure array):
%
%                  Var.name  => variable name (string)
%
%                  Var.type  => external data type (string or numeric):
%                               'byte'   or  nc_byte
%                               'char'   or  nc_char
%                               'short'  or  nc_short
%                               'int'    or  nc_int
%                               'float'  or  nc_float
%                               'double' or  nc_double
%                                     
%                  Var.dimid => dimension IDs (numeric), if not present
%                               or empthy [] the variable is a scalar
%                  
%                  Any other field in the structure, if any,
%                  is proccessed as a variable attribute. For
%                  example:
%
%                  Var.long_name  => "long_name" attribute (string)
%                  Var.add_offset => "add_offset" attribute (number)
%
%                  the values of the attribute can be numeric (scalar
%                  or vector) or characters (array or cell array).
%
% On Output:
%
%    varid       Variable ID
%
%    status      Error flag
%
  
% Initialize.

status = 0;

% Check input structure.

if (~isfield(Var,'name')),
  disp(' ');
  error('NC_VDEF_MATLAB - Cannot field ''name'' in structure array: Var');
end

if (~isfield(Var,'type')),
  disp(' ');
  error('NC_VDEF_MATLAB - Cannot field ''type'' in structure array: Var');
end

if (~isfield(Var,'dimid')),
  disp(' ');
  error('NC_VDEF_MATLAB - Cannot field ''dimid'' in structure array: Var');
end

%--------------------------------------------------------------------------
% Get some NetCDF parameters.
%--------------------------------------------------------------------------

% Set variable external data type representation.

if (isfield(Var,'type')),
  type = Var.type;
  if (ischar(type)),
    switch type
      case 'byte'
        vartyp = netcdf.getConstant('nc_byte');
      case 'char'
        vartyp = netcdf.getConstant('nc_char');
      case 'short'
        vartyp = netcdf.getConstant('nc_short');
      case 'int'
        vartyp = netcdf.getConstant('nc_int');
      case 'float'
        vartyp = netcdf.getConstant('nc_float');
      case 'double'
        vartyp = netcdf.getConstant('nc_double');
    end
  else
    vartyp = type;
  end
else
  error(['NC_VDEF_MATLAB: external data type field ''type'' is ',       ...
	 'missing in structure: Var']);
end

%--------------------------------------------------------------------------
% Define requested variable. Notice that the dimensions ID are flipped to
% column-major order as in Fortran.
%--------------------------------------------------------------------------

if (isfield(Var,'name') && isfield(Var,'dimid')),
  vdid = fliplr(Var.dimid);
  varid = netcdf.defVar(ncid, Var.name, vartyp, vdid);
end

%--------------------------------------------------------------------------
% Add variable attributes.
%--------------------------------------------------------------------------

names = fieldnames(Var);
nfields = length(names);

for n = 1:nfields,
  Aname = char(names(n));
  switch Aname
    case {'name', 'type', 'dimid'}
      put_attribute = false;
    case {'FillValue'}
      put_attribute = true;
      Aname = '_FillValue';
    otherwise,
      put_attribute = true;
      value = Var.(Aname);    
  end

% Define variable attributes.  

  if (put_attribute),

% Attribute value is character array or cell array.

    if (iscellstr(value) || ischar(value)),
      if (iscellstr(value)),
        value = char(value);
      end
      netcdf.putAtt(ncid, varid, Aname, value);

% Attribute value is numeric (scalar or vector). Check external data
% representation.

    else

      switch (vartyp)
        case (netcdf.getConstant('nc_byte'))
          value = int8(value);
        case (netcdf.getConstant('nc_short'))
          value = int16(value);
        case (netcdf.getConstant('nc_int'))   
          value = int32(value);
        case (netcdf.getConstant('nc_float'))
          value = single(value);
        case (netcdf.getConstant('nc_double'))
          value = double(value);
      end
      netcdf.putAtt(ncid, varid, Aname, value);    
    end
  end
end

return

%--------------------------------------------------------------------------

function [varid, status] = nc_vdef_mexnc(ncid, Var)

%
% NC_VDEF_MEXNC:   Create a ROMS variable in a NetCDF file
%
% [varid, status] = nc_vdef(ncid, Var)
%
% This function defines a variable into NetCDF file using the MEXNC
% inteface.
%
% On Input:
%
%    ncid        NetCDF file ID
%
%    Var         Variable information (structure array):
%
%                  Var.name  => variable name (string)
%
%                  Var.type  => external data type (string or numeric):
%                               'byte'   or  nc_byte
%                               'char'   or  nc_char
%                               'short'  or  nc_short
%                               'int'    or  nc_int
%                               'float'  or  nc_float
%                               'double' or  nc_double
%                                     
%                  Var.dimid => dimension IDs (numeric), if not present
%                               or empthy [] the variable is a scalar
%                  
%                  Any other field in the structure, if any,
%                  is proccessed as a variable attribute. For
%                  example:
%
%                  Var.long_name  => "long_name" attribute (string)
%                  Var.FillValue  => "_FillValue" attribute (number)
%
%                  the values of the attribute can be numeric (scalar
%                  or vector) or characters (array or cell array).
%
% On Output:
%
%    varid       Variable ID
%
%    status      Error flag
%
  
% Check input structure.

if (~isfield(Var,'name')),
  disp(' ');
  error('NC_VDEF_MEXNC - Cannot field ''name'' in structure array: Var');
end

if (~isfield(Var,'type')),
  disp(' ');
  error('NC_VDEF_MEXNC - Cannot field ''type'' in structure array: Var');
end

if (~isfield(Var,'dimid')),
  disp(' ');
  error('NC_VDEF_MEXNC - Cannot field ''dimid'' in structure array: Var');
end

%--------------------------------------------------------------------------
% Get some NetCDF parameters.
%--------------------------------------------------------------------------

% Set variable external data type representation.

if (isfield(Var,'type')),
  type = Var.type;
  if (ischar(type)),
    switch type
      case 'byte'
        vartyp = nc_constant('nc_byte');
      case 'char'
        vartyp = nc_constant('nc_char');
      case 'short'
        vartyp = nc_constant('nc_short');
      case 'int'
        vartyp = nc_constant('nc_int');
      case 'float'
        vartyp = nc_constant('nc_float');
      case 'double'
        vartyp = nc_constant('nc_double');
    end
  else
    vartyp = type;
  end
else
  error(['NC_VDEF_MEXNC: external data type field ''type'' is ',        ...
	 'missing in structure: Var']);
end

%--------------------------------------------------------------------------
% Define requested variable.
%--------------------------------------------------------------------------

if (isfield(Var,'name') && isfield(Var,'dimid')),
  vdid = Var.dimid;
  nvdim = length(vdid);
  [varid,status] = mexnc('def_var',ncid,Var.name,vartyp,nvdim,vdid);
  if (varid == -1 || status ~= 0),
    disp(' ');
    disp(mexnc('strerror',status));
    error(['NC_VDEF_MEXNC: DEF_VAR - unable to define variable: ',      ...
           Var.name]);
  end
end

%--------------------------------------------------------------------------
% Add variable attributes.
%--------------------------------------------------------------------------

names = fieldnames(Var);
nfields = length(names);

for n = 1:nfields,
  Aname = char(names(n));
  switch Aname
    case {'name', 'type', 'dimid'}
      put_attribute = false;
    case {'FillValue'}
      put_attribute = true;
      Aname = '_FillValue';
    otherwise,
      put_attribute = true;
      value = Var.(Aname);    
  end

% Define variable attributes.  

  if (put_attribute),

% Attribute value is character array or cell array.

    if (iscellstr(value) || ischar(value)),
      if (iscellstr(value)),
        value = char(value);          % need figure out this one
      end,
      lstr = length(value);
      status = mexnc('put_att_text',ncid,varid,Aname,                   ...
                     nc_constant('nc_char'),lstr,value);
      if (status ~= 0),
        disp(' ');
        disp(mexnc('strerror',status));
        error(['NC_VDEF_MEXNC: PUT_ATT_TEXT - unable to define ',       ...
               'attribute: ',Var.name,':',Aname,'.']);
      end

% Attribute value is numeric (scalar or vector). Check external data
% representation.

    else

      nval = length(value);
      switch (vartyp)
        case (nc_constant('nc_int'))   
          value = int32(value);
          status = mexnc('put_att_int',   ncid,varid,Aname,vartyp,      ...
                         nval,value);
          if (status ~= 0),
            disp(' ');
            disp(mexnc('strerror',status));
            error(['NC_VDEF_MEXNC: PUT_ATT_INT - unable to define ',    ...
                   'attribute: ',Vname.name,':',Aname,'.']);
          end
        case (nc_constant('nc_float'))
          value = single(value);
          status = mexnc('put_att_float', ncid,varid,Aname,vartyp,      ...
                          nval,value);
          if (status ~= 0),
            disp(' ');
            disp(mexnc('strerror',status));
            error(['NC_VDEF_MEXNC: PUT_ATT_FLOAT - unable to define ',  ...
                   'attribute: ',Vname.name,':',Aname,'.']);
          end
        case (nc_constant('nc_double'))
          value = double(value);
          status = mexnc('put_att_double',ncid,varid,Aname,vartyp,      ...
                         nval,value);
          if (status ~= 0),
            disp(' ');
            disp(mexnc('strerror',status));
            error(['NC_VDEF_MEXNC: PUT_ATT_DOUBLE - unable to define ', ...
                   'attribute: ',Vname.name,':',Aname,'.']);
          end
      end
    end
  end
end

return
