function status = nc_attadd(ncfile, Aname, Avalue, varargin)

%
% NC_ATTADD:  Add/modify a global or variable NetCDF attribute
%
% status = nc_attadd(ncfile, Aname, Avalue, Vname)
%
% This function adds/modify a global or variable attribute in a NetCDF
% file. If the "Vname" argument is missing, it is assumed that "Aname"
% is a global attribute.
%
% On Input:
%
%    ncfile     NetCDF file name (string)
%
%    Aname      Attribute name (string)
%
%    Avalue     Attribute value (numeric or string)
%
%    Vname      Variable name (string; optional)
%
% On Output:
%
%    status     Error flag
%

% svn $Id: nc_attadd.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

Vname = [];
switch numel(varargin)
  case 1
    Vname = varargin{1};
end

% Inquire about the contents of the NetCDF file.

Info = nc_inq(ncfile);

% Choose NetCDF file interface.

[method,~,~] = nc_interface(ncfile);

switch(method),
  case {'native'}
    status = nc_attadd_matlab(ncfile, Aname, Avalue, Vname, Info);
  case {'mexnc'}
    status = nc_attadd_mexnc (ncfile, Aname, Avalue, Vname, Info);
  case {'java'}
    error('NC_ATTADD: it is not possible to write into an OpenDAP file.');
  otherwise
    error('NC_ATTADD: unable to determine NetCDF processing interface.');
end

return

%--------------------------------------------------------------------------

function status = nc_attadd_matlab(ncfile, Aname, Avalue, Vname, Info)

%
% NC_ATTADD_MEXNC:  Add/modify a global or variable NetCDF attribute
%
% status = nc_attadd_mexnc(ncfile, Aname, Avalue, Vname, Info)
%
% This function adds/modify a global or variable in a NetCDF file using
% the native Matlab interface. If the "Vname" argument is empty, it is
% assumed that "Aname" is a global attribute.
%
% On Input:
%
%    ncfile      NetCDF file name (string)
%
%    Aname      Attribute name (string)
%
%    Avalue     Attribute value (numeric or string)
%
%    Vname      Variable name (string)
%
%    Info       NetCDF information structure (struct array)
%
% On Output:
%
%    status     Error flag
%

% Initialize error flag.

status = 0;

if (isempty(Vname)),
  got_var = false;
else
  got_var = true;
end

% Open NetCDF file and put it into define mode.

ncid  = netcdf.open(ncfile, 'nc_write');
netcdf.reDef(ncid);

%---------------------------------------------------------------------------
% Add/modify a variable attribute.
%---------------------------------------------------------------------------

if (got_var),

% Get variable ID.

  if (~any(strcmp({Info.Variables.Name}, Vname)))
    nc_inq(ncfile, true);
    disp(' ');
    error(['NC_ATTADD_MEXNC: cannot find NetCDF variable: ',Vname]);
  else
    varid = netcdf.inqVarID(ncid, Vname);
    ivar  = strcmp({Info.Variables.Name}, Vname);
    vtype = Info.Variables(ivar).ncType;  
  end

% Inquire if variable attribute exist.

  if (any(strcmp({Info.Variables(ivar).Attributes.Name}, Aname)))
    disp(' ');
    disp(['Requested attribute "',Aname,'" already exist in ',          ...
          'variable "',Vname,'".']);
    if ischar(Avalue),
      disp(['Updating its value to: "',Avalue,'".']);
    else
      disp(['Updating its value to: ',num2str(Avalue),'.']);
    end
  end

% Add/modify variable attribute.

  if ischar(Avalue),
    netcdf.putAtt(ncid, varid, Aname, Avalue);
  else
    if (strcmpi(Aname,'_FillValue')    ||                               ...
        strcmpi(Aname,'missing_value'))
      switch (vtype)
        case (netcdf.getConstant('nc_byte'))
          value = int8(Avalue);
        case (netcdf.getConstant('nc_short'))
          value = int16(Avalue);
        case (netcdf.getConstant('nc_int'))
          value = int32(Avalue);
        case (netcdf.getConstant('nc_float'))
          value = single(Avalue);
        case (netcdf.getConstant('nc_double'))
          value = double(Avalue);
      end
      netcdf.putAtt(ncid, varid, Aname, value);
    else
      netcdf.putAtt(ncid, varid, Aname, Avalue);
    end    
  end

%---------------------------------------------------------------------------
%  Add/modify a global attribute.
%---------------------------------------------------------------------------
   
else

% Inquire if variable attribute exist.

  if (any(strcmp({Info.Attributes.Name}, Aname)))
    disp(' ');
    disp(['Requested attribute "',Aname,'" already exist in file: ',    ...
          ncfile,'.']);
    if ischar(Avalue),
      disp(['Updating its value to: "',Avalue,'".']);
    else
      disp(['Updating its value to: ',num2str(Avalue),'.']);
    end
  end
  
%  Add/modify character attribute.

  varid  = netcdf.getConstant('nc_global');
  netcdf.putAtt(ncid, varid, Aname, Avalue);
  
end

% Exit definition mode and close NetCDF file.

netcdf.endDef(ncid);
netcdf.close(ncid);

return

%--------------------------------------------------------------------------

function status = nc_attadd_mexnc(ncfile, Aname, Avalue, Vname, Info)

%
% NC_ATTADD_MEXNC:  Add/modify a global or variable NetCDF attribute
%
% status = nc_attadd_mexnc(ncfile, Aname, Avalue, Vname, Info)
%
% This function adds/modify a global or variable in a NetCDF file using
% the MEXNC interface. If the "Vname" argument is empty, it is assumed
% that "Aname" is a global attribute.
%
% On Input:
%
%    ncfile      NetCDF file name (string)
%
%    Aname      Attribute name (string)
%
%    Avalue     Attribute value (numeric or string)
%
%    Vname      Variable name (string)
%
%    Info       NetCDF information structure (struct array)
%
% On Output:
%
%    status     Error flag
%

% Initialize error flag.

if (isempty(Vname)),
  got_var = false;
else
  got_var = true;
end

% Open NetCDF file.

ncid = mexnc('open',ncfile,'nc_write');
if (ncid < 0),
  disp(' ');
  error(['NC_ATTADD_MEXNC: open - unable to open file: ', ncfile]);
end

% Put open file into define mode.

status = mexnc('redef',ncid);
if (status < 0),
  status = -1;
  disp(' ');
  disp(mexnc('strerror',status));
  mexnc('close',ncid);
  error(['NC_ATTADD_MEXNC: redef - unable to put in definition mode: ',  ...
         ncfile]);
end

%---------------------------------------------------------------------------
% Add/modify a variable attribute.
%---------------------------------------------------------------------------

if (got_var),

% Get variable ID.

  if (~any(strcmp({Info.Variables.Name}, Vname)))
    nc_inq(ncfile, true);
    disp(' ');
    error(['NC_ATTADD_MEXNC: cannot find NetCDF variable: ',Vname]);
  end
  varid = mexnc('inq_varid',ncid,Vname);

% Inquire number of variable attributes.

  [nvatts,status] = mexnc('inq_varnatts',ncid,varid);
  if (status < 0),
    disp(' ');
    disp(mexnc('strerror',status));
    error(['NC_ATTADD_MEXNC: inq_varnatts - unable to inquire number ', ...
           'of variable attributes: ',Vname]);
  end

% Inquire if variable attribute exist.

  found = false;
  
  for i=0:nvatts-1
    [attnam,status] = mexnc('inq_attname',ncid,varid,i);
    if (status < 0),
      disp(' ');
      disp(mexnc('strerror',status));
      error(['NC_ATTADD_MEXNC: inq_attname: error while inquiring ',    ...
             'attribute: ',num2str(i)]);
    end,
    if (strcmp(Aname, attnam)),
      found = true;
      break
    end
  end
  if (found),
    disp(' ');
    disp(['Requested attribute "',Aname,'" already exist in ',          ...
          'variable "',Vname,'".']);
    if ischar(Avalue),
      disp(['Updating its value to: "',Avalue,'".']);
    else
      disp(['Updating its value to: ',num2str(Avalue),'.']);
    end
  end

% Inquire variable type.

  [vtype,status] = mexnc('inq_vartype',ncid,varid);
  if (status < 0),
    disp(' ');
    disp(mexnc('strerror',status));
    error(['NC_ATTADD_MEXNC: inq_vartype - unable to inquire ',         ...
           'datatype for variable: ',Vname]);
  end

% Add/modify variable attribute.

  if ischar(Avalue),
    lstr = length(Avalue);
    status = mexnc('put_att_text',ncid,varid,Aname,                     ...
                   nc_constant('nc_char'),lstr,Avalue);
    if (status < 0),
      disp(' ');
      disp(mexnc('strerror',status));
      error(['NC_ATTADD_MEXNC: put_att_text - unable to define ',       ...
             'attribute "',Aname,'"in variable: ',Vname,'.']);
    end
  else
    nval = length(Avalue);
    switch (vtype)
      case (nc_constant('nc_int'))
        value = int32(Avalue);
        status = mexnc('put_att_int',   ncid,varid,Aname,vtype,nval,value);
        if (status < 0),
          disp(' ');
          disp(mexnc('strerror',status));
          error(['NC_ATTADD_MEXNC: put_att_int - unable to define ',    ...
                 'attribute "',Aname,'"in variable: ',Vname,'.']);
        end
      case (nc_constant('nc_float'))
        value = single(Avalue);
        status = mexnc('put_att_float', ncid,varid,Aname,vtype,nval,value);
        if (status < 0),
          disp(' ');
          disp(mexnc('strerror',status));
          error(['NC_ATTADD_MEXNC: put_att_float - unable to define ',  ...
                 'attribute "',Aname,'"in variable: ',Vname,'.']);
        end
      case (nc_constant('nc_double'))
        value = double(Avalue);
        status = mexnc('put_att_double',ncid,varid,Aname,vtype,nval,value);
        if (status < 0),
          disp(' ');
          disp(mexnc('strerror',status));
          error(['NC_ATTADD_MEXNC: put_att_double - unable to define ', ...
                 'attribute "',Aname,'"in variable: ',Vname,'.']);
        end
    end
  end

%---------------------------------------------------------------------------
%  Add/modify a global attribute.
%---------------------------------------------------------------------------
   
else

%  Inquire number of global attributes.

  [natts,status]=mexnc('inq_natts',ncid);
  if (status < 0),
    disp(' ');
    disp(mexnc('strerror',status));
    error(['NC_ATTADD_MEXNC: inq_natts - unable to inquire number of ', ...
           'global attributes: ',ncfile]);
  end
  
%  Inquire if requested global attribute exist.

  found = false;
  
  for i=0:natts-1
    [attnam,status] = mexnc('inq_attname',ncid,                         ...
                            nc_constant('nc_global'),i);
    if (status < 0),
      disp(' ');
      disp(mexnc('strerror',status));
      error(['NC_ATTADD_MEXNC: inq_attname: error while inquiring ',    ...
             'attribute: ',num2str(i)]);
    end
    if (strcmp(Aname, attnam)),
      found = true;
      break
    end
  end
  if (found),
    disp(' ');
    disp(['Requested attribute "',Aname,'" already exist in file: ',    ...
          ncfile,'.']);
    disp(['Updating its value to: "',Avalue,'".']);
  end
  
%  Add/modify character attribute.

  if (ischar(Avalue)),
    lstr = length(Avalue);
    status = mexnc('put_att_text',ncid,nc_constant('nc_global'),Aname,  ...
                   nc_constant('nc_char'),lstr,Avalue);
    if (status < 0),
      disp(' ');
      disp(mexnc('strerror',status));
      error(['NC_ATTADD_MEXNC: put_att_text - unable to define ',       ...
             'attribute "',Aname,'"in file: ',ncfile,'.']);
    end,
  elseif (isinteger(Avalue)),
    nval = length(Avalue); 
    status = mexnc('put_att_int',ncid,nc_constant('nc_global'),Aname,   ...
                   nc_constant('nc_int'),nval,Avalue);
    if (status ~= 0),
      disp(' ');
      disp(mexnc('strerror',status));
      error(['NC_ATTADD_MEXNC: PUT_ATT_INT - unable to define ',        ...
             'attribute: "',Aname,'"in file: ',ncfile,'.']);
    end
  elseif (isfloat(Avalue)),
    nval = length(Avalue);
    status = mexnc('put_att_double',ncid,nc_constant('nc_global'),Aname,...
                   nc_constant('nc_double'),nval,double(Avalue));
    if (status ~= 0),
      disp(' ');
      disp(mexnc('strerror',status));
      error(['NC_ATTADD_MEXNC: PUT_ATT_FLOAT - unable to define ',      ...
             'attribute: "',Aname,'"in file: ',ncfile,'.']);
    end
  end
  
end

% Exit definition mode.

status = mexnc('enddef',ncid);
if (status < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  error(['NC_ATTADD_MEXNC: enddef - unable to exit definition mode: ',  ...
         ncfile]);
end,

%  Close NetCDF file.

cstatus = mexnc('ncclose',ncid);
if (cstatus < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  error(['NC_ATTADD_MEXNC: ncclose - unable to close NetCDF file: ',    ...
         ncfile]);
end

return
