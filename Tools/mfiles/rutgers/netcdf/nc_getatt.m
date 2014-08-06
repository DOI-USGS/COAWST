function A = nc_getatt(ncfile, Aname, Vname)

%
% NC_GETATT:  Gets a global or variable NetCDF attribute
%
% A = nc_getatt(ncfile, Aname,Vname)
%
% This function reads the global or variable attributes from input
% NetCDF file. If the "Vname" argument is missing, it is assumed
% that "Aname" is a global attribute. If both "Aname" and "Vname"
% arguments are missing, it return all global attributes.
%
% On Input:
%
%    ncfile     NetCDF file name or URL file name (string)
%    Aname      Attribute name (string; optional)
%    Vname      Variable name  (string; optional)
%
% On Output:
%
%    A          Attribute value (numeric, string, or struct array)
%
% Usage:
%
%    A = nc_getatt ('ocean_his.nc')
%
%        will return a structure array with all the global attributes
%        names and values: A.Name(:), A.Value(:).
%
%    A = nc_getatt ('ocean_his.nc', 'history')
%
%        will return a string with the value of the global history
%        attribute.
%
%    A = nc_getatt ('ocean_his.nc', [], 'temp')
%
%        will return a structure array with all the 'temp' variable
%        attributes names and values: A.Name(:), A.Value(:).
%
%    A = nc_getatt ('ocean_his.nc', 'FillValue', 'temp')
%
%        will return a numeric value for the attribute 'FillValue'
%        in variable 'temp'.
  
% svn $Id: nc_getatt.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

if (nargin < 3),
  Vname=[];
end
if (nargin < 2),
  Aname=[];
end

% Choose NetCDF file interface.

[method,~,~] = nc_interface(ncfile);

switch(method),
  case {'native'}
    A = nc_getatt_matlab(ncfile, Aname, Vname);  % Matlab native interface
  case {'java'}
    A = nc_getatt_java  (ncfile, Aname, Vname);  % SNCTOOLS JAVA interface
  case {'mexnc'}
    A = nc_getatt_mexnc (ncfile, Aname, Vname);  % MEXNC inteface
  otherwise
    error('NC_GETATT: unable to determine NetCDF processing interface');
end

return

%--------------------------------------------------------------------------

function A = nc_getatt_matlab(ncfile, Aname, Vname)

%
% NC_GETATT_MATLAB:  Gets a global or variable NetCDF attribute
%
% A = nc_getatt_matlab(ncfile, Aname, Vname)
%
% This function reads the global or variable attributes from input
% NetCDF file. It uses the native matlab function "ncinfo" that it
% is available in version 2012a or higher.
% 

%  Initialize.

if (~isempty(Vname)),
  get_varatt = true;          % gets variable attribute(s)
  get_global = false;
else
  get_varatt = false;
  get_global = true;          % gets global attribute(s)
end
if (isempty(Aname)),
  get_allatt = true;          % gets all attributes
else
  get_allatt = false;         % gets specified attribute
end

if (get_allatt),
  A = struct('Name', '', 'Value', []);
else
  A = [];
end

%  Inquire information from URL NetCDF file.

if (get_varatt),
  Info = ncinfo(ncfile,Vname);
else
  Info = ncinfo(ncfile);
end

%--------------------------------------------------------------------------
% Get variable attribute.
%--------------------------------------------------------------------------

if (get_varatt),
  nvatts = length(Info.Attributes);
  for n=1:nvatts,
    if (get_allatt),
      A(n).Name  = Info.Attributes(n).Name;
      A(n).Value = Info.Attributes(n).Value;
    else
      attname = Info.Attributes(n).Name;
      if (strcmp(attname, Aname)),
        A = Info.Attributes(n).Value;
      end
    end
  end

  if (~get_allatt && isempty(A)),
    disp(' ')
    disp(['   Variable attribute: ''',Aname,                            ...
          ''' not found in ''', Vname,''' ...'])
    disp(' ')
  end
end

%--------------------------------------------------------------------------
% Get global attribute.
%--------------------------------------------------------------------------

if (get_global),
  natts = length(Info.Attributes);
  for n=1:natts,
    if (get_allatt),
      A(n).Name  = Info.Attributes(n).Name;
      A(n).Value = Info.Attributes(n).Value;
    else
      attname = Info.Attributes(n).Name;
      if (strcmp(attname, Aname)),
        A = Info.Attributes(n).Value;
      end
    end
  end
  
  if (~get_allatt && isempty(A)),
    disp(' ')
    disp(['   Global attribute: ''',Aname,''' not found ...'])
    disp(' ')
  end
  
end

return

%--------------------------------------------------------------------------

function A = nc_getatt_java(ncfile, Aname, Vname)

%
% NC_GETATT_JAVA:  Gets a global or variable NetCDF attribute
%
% A = nc_getatt_java(ncfile, Aname, Vname)
%
% This function reads the global or variable attributes from input
% NetCDF file. It uses SNCTOOLS function "nc_info".

%  Initialize.

if (~isempty(Vname)),
  get_varatt = true;          % gets variable attribute(s)
  get_global = false;
else
  get_varatt = false;
  get_global = true;          % gets global attribute(s)
end
if (isempty(Aname)),
  get_allatt = true;          % gets all attributes
else
  get_allatt = false;         % gets specified attribute
end

if (get_allatt),
  A = struct('Name', '', 'Value', []);
else
  A = [];
end

%  Inquire information from URL NetCDF file.

Info = nc_info(ncfile); 

%--------------------------------------------------------------------------
% Get variable attribute.
%--------------------------------------------------------------------------

if (get_varatt),
  nvars = length(Info.Dataset);
  for n=1:nvars,
    varname = Info.Dataset(n).Name;
    if (strcmp(varname, Vname)),
      nvatts = length(Info.Dataset(n).Attribute);
      for m=1:nvatts,
        if (get_allatt),
          A(m).Name  = Info.Dataset(n).Attribute(m).Name;
          A(m).Value = Info.Dataset(n).Attribute(m).Value;
        else
          attname = Info.Dataset(n).Attribute(m).Name;
          if (strcmp(attname, Aname)),
            A = Info.Dataset(n).Attribute(m).Value;
          end
        end
      end
    end
  end

  if (~get_allatt && isempty(A)),
    disp(' ')
    disp(['   Variable attribute: ''',Aname,                            ...
          ''' not found in ''', Vname,''' ...'])
    disp(' ')
  end
end

%--------------------------------------------------------------------------
% Get global attribute.
%--------------------------------------------------------------------------

if (get_global),
  natts = length(Info.Attribute);
  for n=1:natts,
    if (get_allatt),
      A(n).Name  = Info.Attribute(n).Name;
      A(n).Value = Info.Attribute(n).Value;
    else
      attname = Info.Attribute(n).Name;
      if (strcmp(attname, Aname)),
        A = Info.Attribute(n).Value;
      end
    end
  end
  
  if (~get_allatt && isempty(A)),
    disp(' ')
    disp(['   Global attribute: ''',Aname,''' not found ...'])
    disp(' ')
  end
  
end

return

%--------------------------------------------------------------------------

function A = nc_getatt_mexnc(ncfile, Aname, Vname)

%
% NC_GETATT_MEXNC:  Gets a global or variable NetCDF attribute
%
% A = nc_getatt_mexnc(ncfile, Aname, Vname)
%
% This function reads the global or variable attributes from input
% NetCDF file. It uses the MEXNC interface.
%

%  Initialize.

if (~isempty(Vname)),
  get_varatt = true;          % gets variable attribute(s)
  get_global = false;
else
  get_varatt = false;
  get_global = true;          % gets global attribute(s)
end
if (isempty(Aname)),
  get_allatt = true;          % gets all attributes
else
  get_allatt = false;         % gets specified attribute
end

if (get_allatt),
  A = struct('Name', '', 'Value', []);
else
  A = [];
end

%  Open NetCDF file.

[ncid]=mexnc('open',ncfile,'nc_nowrite');
if (ncid < 0),
  disp('  ');
  error(['NC_GETATT: open - unable to open file: ', ncfile]);
end

%--------------------------------------------------------------------------
% Get variable attribute(s).
%--------------------------------------------------------------------------

if (get_varatt),

%  Get variable ID.

  [varid]=mexnc('inq_varid',ncid,Vname);
  if (varid < 0),
    [~]=mexnc('close',ncid);
    disp('  ');
    error(['NC_GETATT: inq_varid - cannot find variable: ',Vname]);
  end

%  Inquire number of variable attributes.

  [nvatts,status]=mexnc('inq_varnatts',ncid,varid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_GETATT: inq_varnatts - unable to inquire number of ',    ...
           'variable attributes: ',Vname]);
  end

%  Inquire if variable attribute exist.

  found = false;
  
  for i=1:nvatts
    [attnam,status]=mexnc('inq_attname',ncid,varid,i-1);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_GETATT: inq_attname: error while inquiring ',          ...
             'attribute:', num2str(i-1)]);
    end
    if (get_allatt),
      A(i).Name = attnam;
      found = true;
    else
      if (strcmp(attnam, Aname)),
        found = true;
      end
    end

%  Read in requested attribute.
  
    if (found),
      [atype,status]=mexnc('inq_atttype',ncid,varid,attnam);
      if (status < 0),
        disp('  ');
        disp(mexnc('strerror',status));
        error(['NC_GETATT: inq_atttype - unable to inquire datatype ',  ...
               'for attribute: ',attnam]);
      end

      switch (atype)
        case (nc_constant('nc_char'))
          [avalue,status]=mexnc('get_att_text',ncid,varid,attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_GETATT: get_att_text - unable to read ',         ...
                   'attribute "',attnam,'"in variable: ',Vname,'.']);
          end
        case (nc_constant('nc_int'))
          [avalue,status]=mexnc('get_att_int', ncid,varid,attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_GETATT: get_att_int - unable to read ',          ...
                   'attribute "',attnam,'"in variable: ',Vname,'.']);
          end
        case (nc_constant('nc_float'))
          [avalue,status]=mexnc('get_att_float',ncid,varid,attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_GETATT: get_att_float - unable to read ',        ...
                   'attribute "',attnam,'"in variable: ',Vname,'.']);
          end
        case (nc_constant('nc_double'))
          [avalue,status]=mexnc('get_att_double',ncid,varid,attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_GETATT: get_att_double - unable to read ',       ...
                   'attribute "',attnam,'"in variable: ',Vname,'.']);
          end
      end
      if (get_allatt),
	A(i).Value = avalue;
      else
	A = avalue;
	break
      end
    end
  end
  
  if (~get_allatt && isempty(A)),
    disp(' ')
    disp(['   Variable attribute: ''',Aname,     ...
          ''' not found in ''', Vname,''' ...'])
    disp(' ')
  end
end

%--------------------------------------------------------------------------
% Get global attribute(s).
%--------------------------------------------------------------------------

if (get_global),

%  Inquire number of global attributes.

  [natts,status]=mexnc('inq_natts',ncid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_GETATT: inq_natts - unable to inquire number of global', ...
           ' attributes: ',ncfile]);
  end
  
%  Inquire if requested global attribute exist.

  found = false;
  
  for i=1:natts,
    [attnam,status]=mexnc('inq_attname',ncid,nc_constant('nc_global'),  ...
                          i-1);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_GETATT: inq_attname: error while inquiring ',          ...
             'attribute: ', num2str(i-1)]);
    end
    if (get_allatt),
      A(i).Name = attnam;
      found = true;
    else
      if (strcmp(attnam, Aname)),
        found = true;
      end
    end
  
%  Read in requested attribute.

    if (found),
      [atype,status]=mexnc('inq_atttype',ncid,nc_constant('nc_global'), ...
                           attnam);
      if (status < 0),
        disp('  ');
        disp(mexnc('strerror',status));
        error(['NC_GETATT: inq_atttype - unable to inquire datatype ',  ...
               'for attribute: ',attnam]);
      end

      switch (atype)
        case (nc_constant('nc_char'))
          [avalue,status]=mexnc('get_att_text',ncid,                    ...
                                nc_constant('nc_global'),attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_GETATT: get_att_text - unable to read global ',  ...
                   'attribute: "',attnam,'".']);
          end
        case (nc_constant('nc_int'))
          [avalue,status]=mexnc('get_att_int', ncid,                    ...
                                nc_constant('nc_global'),attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_GETATT: get_att_int - unable to read global ',   ...
                   'attribute: "',attnam,'".']);
          end
        case (nc_constant('nc_float'))
          [avalue,status]=mexnc('get_att_float',ncid,                   ...
                                nc_constant('nc_global'),attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_GETATT: get_att_float - unable to read global ', ...
                   'attribute: "',attnam,'".']);
          end
        case (nc_constant('nc_double'))
          [avalue,status]=mexnc('get_att_double',ncid,                  ...
                                nc_constant('nc_global'),attnam);
          if (status < 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_GETATT: get_att_double - unable to read global', ...
                   ' attribute: "',attnam,'".']);
          end
      end
      if (get_allatt),
	A(i).Value = avalue;
      else
	A = avalue;
	break
      end
    end
  end
  
  if (~get_allatt && isempty(A)),
    disp(' ')
    disp(['   Global attribute: ''',Aname,''' not found ...'])
    disp(' ')
  end
end

%  Close NetCDF file.

[cstatus]=mexnc('ncclose',ncid);
if (cstatus < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_GETATT: ncclose - unable to close NetCDF file: ', ncfile]);
end

return
