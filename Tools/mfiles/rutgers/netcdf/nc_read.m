function f = nc_read(ncfile, Vname, Tindex, ReplaceValue, PreserveType)

%
% NC_READ:  Read requested NetCDF variable
%
% f = nc_read(ncfile, Vname, Tindex, ReplaceValue, PreserveType)
%
% This function reads in a generic variable from a NetCDF file. If only
% the water points are available, this function fill the land areas with
% zero and returns full fields.
%
% On Input:
%
%    ncfile        NetCDF file name or URL name (string)
%
%    Vname         NetCDF variable name to read (string)
%
%    Tindex        Optional, time record index to read (integer):
%                    If Tindex is provided, only the requested time record
%                    is read when the variable has the unlimited dimension
%                    or the word "time" in any of its dimension names.
%                    Otherwise, provide an empty [] argument.
%
%    ReplaceValue  Optional, value to use when _FillValue or missing_value
%                    attribute is found in variable. If not provided, a
%                    zero value will used.
%                    In some instances, like plotting, it is advantageous
%                    to set ReplaceValue = NaN to better visualize the
%                    land masking or the missing data.
%
%    PreserveType  Switch to preserve numerical data type. If false,
%                    convert numerical data to double precision. It
%                    has no effect on data that is already in double
%                    precision.
%
% On Output:
%
%    f             Field (scalar or array)
%
% calls:           nc_inq
%                  nc_vargetr (SNCTOOLS java interface to OpenDAT)
%                  nc_read_matlab, nc_read_java, nc_read_mexnc (private)
%

% svn $Id: nc_read.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% If the PreserveType switch is not provided, set default to convert data
% to double precision.

if (nargin < 5),
  PreserveType = false;
end

% If ReplaceValue is not provided, use zero as a fill value.

if (nargin < 4),
  ReplaceValue = 0;
end

% If Tindex is not provided, use empty value.

if (nargin < 3)
  Tindex = [];
end

%--------------------------------------------------------------------------
% Read in requested variable.
%--------------------------------------------------------------------------

% Inquire about the contents of the NetCDF file.

Info = nc_inq(ncfile);

if (~any(strcmp({Info.Variables.Name}, Vname)))
  nc_inq(ncfile, true);
  disp(' ');
  error(['NC_READ: cannot find NetCDF variable: ',Vname]);
else
  ivar   = find(strcmp({Info.Variables.Name}, Vname));
  nvdims = length(Info.Variables(ivar).Dimensions);
end


% Choose NetCDF file interface. Use 'native' Matlab interface as default.

[method,~,~] = nc_interface(ncfile);

switch(method),
  case {'native'}
    f = nc_read_matlab(ncfile,Vname,Tindex,ReplaceValue,PreserveType,Info);
  case {'java'}
    f = nc_read_java  (ncfile,Vname,Tindex,ReplaceValue,PreserveType,Info);
  case {'mexnc'}
    f = nc_read_mexnc (ncfile,Vname,Tindex,ReplaceValue,PreserveType,Info);
  otherwise
    error('NC_READ: unable to determine NetCDF processing interface');
end

%--------------------------------------------------------------------------
% If NetCDF file is a ROMS file with just the water data, build the full
% array (land and water points).
%--------------------------------------------------------------------------

if (nvdims > 0),
  if (any(strcmp({Info.Variables(ivar).Dimensions.Name}, 'xy_rho' )) || ...
      any(strcmp({Info.Variables(ivar).Dimensions.Name}, 'xyz_rho')));
    f = land_points (ncfile,Vname,ReplaceValue,Info,f);
  end
end

return

%--------------------------------------------------------------------------

function f = land_points (ncfile, Vname, ReplaceValue, Info, f)

%
% LAND_POINTS:  Adds land points to water data
%
% f = land_points(ncfile, Vname, ReplaceValue, Info, f)
%
% This function adds land points to read water data.  It build the full
% data array with land and water points.
%
% On Input:
%
%    ncfile        NetCDF file name or URL name (string)
%
%    Vname         NetCDF variable name to read (string)
%
%    ReplaceValue  Value to use when _FillValue or missing_value
%                    attribute is found in variable
%
%    Info          NetCDF information structure (struct array)
%
%    f             Water point data (1D array)
%
% On Output:
%
%    f             Full field including land and water points (array)
%

% Choose NetCDF file interface.

[Fmethod,~,~] = nc_interface(ncfile);
  
% Inquire ROMS dimension values.

ndims = length(Info.Dimensions);

for n=1:ndims,
  dname =char(Info.Dimensions(n).Name);
  switch dname
    case 'xi_rho',
      Lr = Info.Dimensions(n).Length;
    case 'xi_u',
      Lu = Info.Dimensions(n).Length;
    case 'xi_v',
      Lv = Info.Dimensions(n).Length;
    case 'eta_rho',
      Mr = Info.Dimensions(n).Length;
    case 'eta_u',
      Mu = Info.Dimensions(n).Length;
    case 'eta_v',
      Mv = Info.Dimensions(n).Length;
    case 's_rho',
      Nr = Info.Dimensions(n).Length;
    case 's_w',
      Nw = Info.Dimensions(n).Length;
  end
end  
 
% Check if data is only available at water points.

is2d = false;
is3d = false;

ivar = find(strcmp({Info.Variables.Name}, Vname));

if (~isempty(ivar)),
  nvdims = length(Info.Variables(ivar).Dimensions);

  for n=1:nvdims,
    dname = char(Info.Variables(ivar).Dimensions(n).Name);
    switch dname
      case 'xy_rho',
        msknam = 'mask_rho';
        is2d = true; Im = Lr; Jm = Mr;
      case 'xy_u',
        msknam = 'mask_u';
        is2d = true; Im = Lu; Jm = Mu;
      case 'xy_v',
        msknam = 'mask_v';
        is2d = true; Im = Lv; Jm = Mv;
      case 'xyz_rho',
        msknam = 'mask_rho';
        is3d = true; Im = Lr; Jm = Mr; Km = Nr;
      case 'xyz_u',
        msknam = 'mask_u';
        is3d = true; Im = Lu; Jm = Mu; Km = Nr;
      case 'xyz_v',    
        msknam='mask_v';
        is3d = true; Im = Lv; Jm = Mv; Km = Nr;
      case 'xyz_w',    
        msknam = 'mask_rho';
        is3d = true; Im = Lr; Jm = Mr; Km = Nw;
    end
  end
end

water = is2d | is3d;

% If water data only, read in the appropriate Land/Sea mask.

if (water),
  got_mask = any(strcmp(vnames,msknam));

  if (got_mask),
    switch (Fmethod)
      case {'native'}
        mask = nc_read_matlab(ncfile,msknam,[],ReplaceValue,false,Info);
      case {'java'}
        mask = nc_read_java  (ncfile,msknam,[],ReplaceValue,false,Info);
      case {'mexnc'}
        mask = nc_read_mexnc (ncfile,msknam,[],ReplaceValue,false,Info);
    end
  else  
    Gname = input('Enter grid NetCDF file: ');
    [Gmethod,~,~] = nc_interface(Gname);
    switch (Gmethod)
      case {'native'}
        mask = nc_read_matlab(Gname,msknam,[],ReplaceValue,false,Info);
      case {'java'}
        mask = nc_read_java  (Gname,msknam,[],ReplaceValue,false,Info);
      case {'mexnc'}
        mask = nc_read_mexnc (Gname,msknam,[],ReplaceValue,false,Info);
    end
  end

end

% Build full array with land and water points.

if (water),

  v = f;
  [Npts,Ntime]=size(v);

  if (is2d),
    f = squeeze(ones([Im,Jm,Ntime])).*ReplaceValue;
    MASK = squeeze(repmat(mask,[1,1,Ntime]));
    ind = MASK > 0;
    f(ind) = v;
  elseif (is3d),    
    f = squeeze(ones([Im,Jm,Km,Ntime])).*ReplaceValue;
    MASK = squeeze(repmat(mask,[1,1,Km,Ntime]));
    ind = MASK > 0;
    f(ind) = v;
  end

end

return

%--------------------------------------------------------------------------

function f = nc_read_matlab(ncfile, Vname, Tindex,                       ...
                            ReplaceValue, PreserveType, Info)

%
% NC_READ_MATLAB:  Internal routine to read requested NetCDF variable
%
% f = nc_read_matlab(ncfile,Vname,Tindex,ReplaceValue,PreserveType,Info)
%
% This function reads in requested variable using Matlab's native NetCDF
% inteface.
%
% On Input:
%
%    ncfile        NetCDF file name or URL name (character string)
%
%    Vname         NetCDF variable name to read (character string)
%
%    Tindex        Time record index to read (integer)
%
%    ReplaceValue  Value to use when "_FillValue" or "missing_value"
%                    attribute is found in variable.
%
%    PreserveType  Switch to preserve numerical data type. If false,
%                    convert numerical data to double precision. It
%                    has not effect data is already in double
%                    precision.
%
%    Info          NetCDF information structure (struct array)
%
% On Output:
%
%    f             Field (scalar or array)
%

% Initialize.

got.add_offset    = false;
got.FillValue     = false;
got.missing_value = false;
got.RecDim        = false;
got.scale_factor  = false;

ivar   = find(strcmp({Info.Variables.Name}, Vname));
nvdims = length(Info.Variables(ivar).Dimensions);
nvatts = length(Info.Variables(ivar).Attributes);
nctype = Info.Variables(ivar).ncType;

% Activate switch for reading specific record.

time_rec = false;
if (~isempty(Tindex)),
  time_rec = true;
end,

% Check if there is an unlimited dimension or a time dimension.

if (nvdims > 0),
  for n=1:nvdims,
    dname = char(Info.Variables(ivar).Dimensions(n).Name);
    if (Info.Variables(ivar).Dimensions(n).Unlimited ||                 ...
        ~isempty(strfind(dname,'time'))),
      got.RecDim = true;
      TimeDimName = dname;
    end
  end
end

% Inquire information about the attributes.

if (nvatts > 0),
  for n=1:nvatts,
    aname = char(Info.Variables(ivar).Attributes(n).Name);
    switch aname
      case 'add_offset'
        offset = Info.Variables(ivar).Attributes(n).Value;
        got.add_offset = true;
      case {'_FillValue', '_fillvalue', 'missing_value'}
        spval = Info.Variables(ivar).Attributes(n).Value;
        if (strcmp(aname, 'missing_value')),
          got.missing_value = true;
        else
          got.FillValue = true;
        end
      case 'scale_factor'
        scale = Info.Variables(ivar).Attributes(n).Value;
        got.scale_factor = true;
    end
  end
end

% Set start and count indices to process.

if (nvdims > 0),
  start = zeros([1 nvdims]);
  count = Inf([1 nvdims]); 

  for n=1:nvdims,
    dname = char(Info.Variables(ivar).Dimensions(n).Name);
    start(n) = 0;
    count(n) = Info.Variables(ivar).Dimensions(n).Length;
    if (time_rec && got.RecDim),
      if (strcmp(dname, TimeDimName)),
        start(n) = Tindex-1;
%       start(n) = Tindex;
        count(n) = 1;
      end
    end
  end
end

%--------------------------------------------------------------------------
% Read in requested variable. Use SNCTOOLS function "nc_vargetr" to read
% variable in its native (raw) data type.
%--------------------------------------------------------------------------

% Open NetCDF file.

ncid  = netcdf.open(ncfile, 'nc_nowrite');
varid = netcdf.inqVarID(ncid, Vname);

% Read in data: scalar or multi-dimensional array.

if (nvdims == 0),
  f = netcdf.getVar(ncid, varid);
else
  f = netcdf.getVar(ncid, varid, start, count);
end

% Close NetCDF file.

netcdf.close(ncid);

%--------------------------------------------------------------------------
% Post-process read data.
%--------------------------------------------------------------------------

% Search for fill value or missing values.

ind = [];

if (got.FillValue || got.missing_value),
  if (iscellstr(f) || ischar(f)),
    ind = find(f == spval);
    f(ind) = spval;
  else
    ind = find(abs(f-spval) < 4*eps(double(spval)));
  end
end

% Scale data and/or add offset value.

if (isnumeric(f)),
  if (got.add_offset),
    if (got.scale_factor),
      switch nctype
       case {netcdf.getConstant('nc_float'),                            ...
             netcdf.getConstant('nc_double')}
          f = f.*scale+offset;
        case {nc_int, nc_short, nc_byte}
          f = double(f).*scale+offset;
      end
    else
      switch nctype
        case {netcdf.getConstant('nc_float'),                           ...
              netcdf.getConstant('nc_double')}
          f = f+offset;
        case {netcdf.getConstant('nc_int'),                             ...
              netcdf.getConstant('nc_short'),                           ...
              netcdf.getConstant('nc_byte')}
          f = double(f)+offset;
      end
    end
  elseif (got.scale_factor);
    switch nctype
      case {netcdf.getConstant('nc_float'),                             ...
            netcdf.getConstant('nc_double')}
        f = f.*scale;
      case {netcdf.getConstant('nc_int'),                               ...
            netcdf.getConstant('nc_short'),                             ...
            netcdf.getConstant('nc_byte')}
        f = double(f).*scale;
    end
  end
end

% Set fill values or missing values with specified ReplaceValue.

if (~isempty(ind) && isnumeric(f)),
  f(ind) = ReplaceValue;
end

%--------------------------------------------------------------------------
% If requested and applicable, convert data to double precision.
%--------------------------------------------------------------------------

if (~PreserveType) && ~(iscellstr(f) || ischar(f)),
  f = double(f);
end

return

%--------------------------------------------------------------------------

function f = nc_read_java(ncfile, Vname, Tindex,                         ...
                          ReplaceValue, PreserveType, Info)

%
% NC_READ_JAVA:  Internal routine to read requested NetCDF variable
%
% f=nc_read_java(ncfile,Vname,Tindex,ReplaceValue,PreserveType,Info)
%
% This function reads in requested variable from a URL NetCDF file
% (OpenDAP).  It uses the SCNTOOLS java interface.  The internal
% parameter PRESERVE_FVD is set to TRUE so the data is always
% processed in column-major order (Fortran, Matlab). The data is not
% transposed by the SCNTOOLS interface.
%
% On Input:
%
%    ncfile         NetCDF file name or URL name (character string)
%
%    Vname         NetCDF variable name to read (character string)
%
%    Tindex        Time record index to read (integer)
%
%    ReplaceValue  Value to use when "_FillValue" or "missing_value"
%                    attribute is found in variable.
%
%    PreserveType  Switch to preserve numerical data type. If false,
%                    convert numerical data to double precision. It
%                    has not effect data is already in double
%                    precision.
%
%    Info          NetCDF information structure (struct array)
%
% On Output:
%
%    f             Field (scalar or array)
%

% Initialize.

got.add_offset    = false;
got.FillValue     = false;
got.missing_value = false;
got.RecDim        = false;
got.scale_factor  = false;

ivar   = find(strcmp({Info.Variables.Name}, Vname));
nvdims = length(Info.Variables(ivar).Dimensions);
nvatts = length(Info.Variables(ivar).Attributes);
nctype = Info.Variables(ivar).ncType;

% Check value of persistent switch to process data in column-major or
% row-major order.

if (ispref('SNCTOOLS','PRESERVE_FVD')),
  saved_preserve_fvd = getpref('SNCTOOLS','PRESERVE_FVD');
else
  saved_preserve_fvd = false;             % default value in SNCTOOLS
end

% Set temporarily persistent switch to process array data in column-major
% order.

setpref('SNCTOOLS','PRESERVE_FVD',true);

% Activate switch for reading specific record.

time_rec = false;
if (~isempty(Tindex)),
  time_rec = true;
end,

% Check if there is an unlimited dimension or a time dimension.

if (nvdims > 0),
  for n=1:nvdims,
    dname = char(Info.Variables(ivar).Dimensions(n).Name);
    if (Info.Variables(ivar).Dimensions(n).Unlimited ||                 ...
        ~isempty(strfind(dname,'time'))),
      got.RecDim = true;
      TimeDimName = dname;
    end
  end
end

% Inquire information about the attributes.

if (nvatts > 0),
  for n=1:nvatts,
    aname = char(Info.Variables(ivar).Attributes(n).Name);
    switch aname
      case 'add_offset'
        offset = Info.Variables(ivar).Attributes(n).Value;
        got.add_offset = true;
      case {'_FillValue', '_fillvalue', 'missing_value'}
        spval = Info.Variables(ivar).Attributes(n).Value;
        if (strcmp(aname, 'missing_value')),
          got.missing_value = true;
        else
          got.FillValue = true;
        end
      case 'scale_factor'
        scale = Info.Variables(ivar).Attributes(n).Value;
        got.scale_factor = true;
    end
  end
end

% Set start and count indices to process.

if (nvdims > 0),
  start = zeros([1 nvdims]);
  count = Inf([1 nvdims]); 

  for n=1:nvdims,
    dname = char(Info.Variables(ivar).Dimensions(n).Name);
    start(n) = 0;
    count(n) = Inf;
    if (time_rec && got.RecDim),
      if (strcmp(dname, TimeDimName)),
        start(n) = Tindex-1;
        count(n) = 1;
      end
    end
  end
end

%--------------------------------------------------------------------------
% Read in requested variable. Use SNCTOOLS function "nc_vargetr" to read
% variable in its native (raw) data type.
%--------------------------------------------------------------------------

% Read in data: scalar or multi-dimensional array.

if (nvdims == 0),
  f = nc_vargetr(ncfile,Vname);
else
  f = nc_vargetr(ncfile,Vname,start,count);
end

% Set persistent switch back to user or SNCTOOLS default value.

setpref('SNCTOOLS','PRESERVE_FVD', saved_preserve_fvd);

%--------------------------------------------------------------------------
% Post-process read data.
%--------------------------------------------------------------------------

% Search for fill value or missing values.

ind = [];

if (got.FillValue || got.missing_value),
  if (iscellstr(f) || ischar(f)),
    ind = find(f == spval);
    f(ind) = spval;
  else
    ind = find(abs(f-spval) < 4*eps(double(spval)));
  end
end

%  Scale data and/or add offset value.

if (isnumeric(f)),
  if (got.add_offset),
    if (got.scale_factor),
      switch nctype
        case {nc_constant('nc_float'),                                  ...
              nc_constant('nc_double')}
          f = f.*scale+offset;
        case {nc_constant('nc_int'),                                    ...
              nc_constant('nc_short'),                                  ...
              nc_constant('nc_byte')}
          f = double(f).*scale+offset;
      end
    else
      switch nctype
        case {nc_constant('nc_float'),                                  ...
              nc_constant('nc_double')}
          f = f+offset;
       case {nc_constant('nc_int'),                                     ...
             nc_constant('nc_short'),                                   ...
             nc_constant('nc_byte')}
          f = double(f)+offset;
      end
    end
  elseif (got.scale_factor);
    switch nctype
      case {nc_constant('nc_float'),                                    ...
            nc_constant('nc_double')}
        f = f.*scale;
      case {nc_constant('nc_int'),                                      ...
            nc_constant('nc_short'),                                    ...
            nc_constant('nc_byte')}
        f = double(f).*scale;
    end
  end
end

%  Set fill values or missing values with specified ReplaceValue.

if (~isempty(ind) && isnumeric(f)),
  f(ind) = ReplaceValue;
end

%--------------------------------------------------------------------------
%  Convert data to double precision.
%--------------------------------------------------------------------------

if (~PreserveType) && ~(iscellstr(f) || ischar(f)),
  f = double(f);
end

return

%--------------------------------------------------------------------------

function f = nc_read_mexnc(ncfile, Vname, Tindex,                        ...
                           ReplaceValue, PreserveType, Info)

%
% NC_READ_MEXNC:  Internal routine to read requested NetCDF variable
%
% f=nc_read_mexnc(ncfile,Vname,Tindex,ReplaceValue,PreserveType,Info)
%
% This function reads in a requested variable form NetCDF file using the
% MEXNC interface.
%
%    ncfile         NetCDF file name or URL name (character string)
%
%    Vname         NetCDF variable name to read (character string)
%
%    Tindex        Time record index to read (integer)
%
%    ReplaceValue  Value to use when "_FillValue" or "missing_value"
%                    attribute is found in variable.
%
%    PreserveType  Switch to preserve numerical data type. If false,
%                    convert numerical data to double precision. It
%                    has not effect data is already in double
%                    precision.
%
%    Info          NetCDF information structure (struct array)
%
% On Output:
%
%    f             Field (scalar or array)
%

%  Set-up printing information switch.

global IPRINT

if (isempty(IPRINT)),
  IPRINT = false;
end

%  Activate switch for reading specific record.

time_rec = false;
if (~isempty(Tindex)),
  time_rec = true;
end,

% Open NetCDF file.

[ncid]=mexnc('ncopen',ncfile,'nc_nowrite');
if (ncid == -1),
  error(['NC_READ_MEXNC: ncopen - unable to open file: ' ncfile])
end

% Supress all error messages from NetCDF.

mexnc('setopts',0);

%--------------------------------------------------------------------------
% Inquire about requested variable.
%--------------------------------------------------------------------------

% Get variable ID.

[varid] = mexnc('ncvarid',ncid,Vname);
if (varid < 0),
  mexnc('ncclose',ncid);
  nc_inq(ncfile);
  disp('  ');
  error(['NC_READ_MEXNC: ncvarid - cannot find variable: ',Vname])
end

% Inquire about unlimited dimension.

[~,~,~,recdim,status] = mexnc('ncinquire',ncid);
if (status == -1),
  error(['NC_READ_MEXNC: ncinquire - cannot inquire file: ',ncfile])
end

% Get information about requested variable.

[Vname,nctype,nvdims,dimids,nvatts,status] = mexnc('ncvarinq',ncid,varid);
if (status == -1),
  error(['NC_READ_MEXNC: ncvarinq - unable to inquire about variable: ',Vname])
end

% Inquire about the _FillValue attribute.

got.add_offset    = false;
got.FillValue     = false;
got.missing_value = false;
got.scale_factor  = false;

for i = 0:nvatts-1,
  [attnam,status] = mexnc('inq_attname',ncid,varid,i);
  if (status == -1)
    error(['NC_READ_MEXNC: inq_attname: error while inquiring ',        ...
           'attribute', num2str(i)]);
  end
  [atype,status] = mexnc('inq_atttype',ncid,varid,attnam);
  if (status == -1)
    error(['NC_READ_MEXNC: inq_atttype: error while inquiring ',        ...
           'attribute ',num2str(i)]);
  end
  switch (attnam)
    case 'add_offset'
      if (atype == nc_constant('nc_double')),
        [offset,status] = mexnc('get_att_double',ncid,varid,attnam); 
      elseif (atype == nc_constant('nc_float')),
        [offset,status] = mexnc('get_att_float' ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_int')),
        [offset,status] = mexnc('get_att_int'   ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_short')),
        [offset,status] = mexnc('get_att_short' ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_byte')),
        [offset,status] = mexnc('get_att_schar' ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_char')),
        [offset,status] = mexnc('get_att_text'  ,ncid,varid,attnam);
      else
        [offset,status] = mexnc('ncattget'      ,ncid,varid,attnam);
      end
      if (status == -1),
        error(['NC_READ_MEXNC: ncattget error while reading: ',         ...
               attnam])
      end
      got.add_offset = true;
    case {'_FillValue', '_fillvalue', 'missing_value'}
      if (atype == nc_constant('nc_double')),
        [spval,status] = mexnc('get_att_double',ncid,varid,attnam); 
      elseif (atype == nc_constant('nc_float')),
        [spval,status] = mexnc('get_att_float' ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_int')),
        [spval,status] = mexnc('get_att_int'   ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_short')),
        [spval,status] = mexnc('get_att_short' ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_byte')),
        [spval,status] = mexnc('get_att_schar' ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_char')),
        [spval,status] = mexnc('get_att_text'  ,ncid,varid,attnam);
      else
        [spval,status] = mexnc('ncattget'      ,ncid,varid,attnam);
      end
      if (status == -1),
        error(['NC_READ_MEXNC: ncattget error while reading: ', attnam]);
      end
      if (strcmp(attnam,'missing_value')),
        got.missing_value = true;
      else
        got.FillValue = true;
      end
    case 'scale_factor'
      if (atype == nc_constant('nc_double')),
        [scale,status] = mexnc('get_att_double',ncid,varid,attnam); 
      elseif (atype == nc_constant('nc_float')),
        [scale,status] = mexnc('get_att_float' ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_int')),
        [scale,status] = mexnc('get_att_int'   ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_short')),
        [scale,status] = mexnc('get_att_short' ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_byte')),
        [scale,status] = mexnc('get_att_schar' ,ncid,varid,attnam);
      elseif (atype == nc_constant('nc_char')),
        [scale,status] = mexnc('get_att_text'  ,ncid,varid,attnam);
      else
        [scale,status] = mexnc('ncattget'      ,ncid,varid,attnam);
      end
      if (status == -1),
        error(['NC_READ_MEXNC: ncattget error while reading: ', attnam]);
      end
      got.scale_factor=true;
  end
end

% Inquire about dimensions.

index = 0;
for n=1:nvdims,
  [name,dsize,status]=mexnc('ncdiminq',ncid,dimids(n));
  if (status == -1),
    error(['NC_READ_MEXNC: ncdiminq - unable to inquire about ',        ...
           'dimension ID: ',num2str(dimids(n))])
  else
    lstr = length(name);
    dimnam(n,1:lstr) = name(1:lstr);
    dimsiz(n) = dsize;
    start(n)  = 0;
    count(n)  = dsize;
    if ((dimids(n) == recdim) || ~isempty(strfind(name,'time'))),
      index = n;
    end
  end
end

%  It reading specific time record, reset variable bounds.

nvdim = nvdims;
if (time_rec && (index > 0)),
  start(index) = Tindex-1;
  count(index) = 1;
  nvdims = nvdims-1;
end

%--------------------------------------------------------------------------
% Read in requested variable.
%--------------------------------------------------------------------------

%  Read in scalar.

if (nvdim == 0),
  switch nctype
    case (nc_constant('nc_double'))
      [f,status] = mexnc('get_var_double',ncid,varid);
      myfunc = 'get_var_double';
    case (nc_constant('nc_float'))
      [f,status] = mexnc('get_var_float' ,ncid,varid);
      myfunc = 'get_var_float';
    case (nc_constant('nc_int'))
      [f,status] = mexnc('get_var_int'   ,ncid,varid);
      myfunc = 'get_var_int';
    case (nc_constant('nc_short'))
      [f,status] = mexnc('get_var_short' ,ncid,varid);
      myfunc = 'get_var_short';
    case (nc_constant('nc_byte'))
      [f,status] = mexnc('get_var_schar' ,ncid,varid);
      myfunc = 'get_var_schar';
    case (nc_constant('nc_char'))
      [f,status] = mexnc('get_var_text'  ,ncid,varid);
      myfunc = 'get_var_text';
    otherwise
      [f,status] = mexnc('ncvarget1'     ,ncid,varid,0);
      myfunc = 'ncvarget1';
  end
  if (status == -1),
    error(['NC_READ_MEXNC: ',myfunc,' - error while reading: ',Vname])
  end

%  Read in a multidemensional array.

else

  switch nctype
    case (nc_constant('nc_double'))
      [f,status] = mexnc('get_vara_double',ncid,varid,start,count);
      myfunc = 'get_vara_double';
    case (nc_constant('nc_float'))
      [f,status] = mexnc('get_vara_float' ,ncid,varid,start,count);
      myfunc = 'get_vara_float';
    case (nc_constant('nc_int'))
      [f,status] = mexnc('get_vara_int'   ,ncid,varid,start,count);
      myfunc = 'get_vara_int';
    case (nc_constant('nc_short'))
      [f,status] = mexnc('get_vara_short' ,ncid,varid,start,count);
      myfunc = 'get_vara_short';
    case (nc_constant('nc_byte'))
      [f,status] = mexnc('get_vara_schar' ,ncid,varid,start,count);
      myfunc = 'get_vara_schar';
    otherwise
      [f,status] = mexnc('ncvarget'       ,ncid,varid,start,count);
      myfunc = 'ncvarget';
  end
  if (status == -1),
    error(['NC_READ_MEXNC: ',myfunc,' - error while reading: ',Vname])
  end

  if (nvdims == 3),
    if (length(start) == 3),
      f = reshape(f,[count(3),count(2),count(1)]);
    elseif (length(start) == 4),
      f = reshape(f,[count(4),count(3),count(2)]);
    end
  end
  
  if (nvdims == 4),
    if (length(start) == 4),
      f = reshape(f,[count(4),count(3),count(2),count(1)]);
    elseif (length(start) == 5),
      f = reshape(f,[count(5),count(4),count(3),count(2)]);
    end  
  end
  
end

% Print information about variable.

if (IPRINT),
  if (nvdims > 0),
    disp(' ')
    disp([Vname ' has the following dimensions (input order):']);
    disp(' ')
    for n=1:nvdim,
      s=[blanks(11) int2str(n) ') ' dimnam(n,:) ' = ' int2str(dimsiz(n))];
      disp(s);
    end
    disp(' ')
    disp([Vname ' loaded into an array of size:  [' int2str(size(f)) ']']);
    disp(' ')
  else
    disp(' ')
    disp([Vname ' is a scalar and has a value of ',num2str(f)]);
    disp(' ')
  end
end

%--------------------------------------------------------------------------
%  Post-process read data.
%--------------------------------------------------------------------------

%  Search for fill value or missing values.

ind = [];

if (got.FillValue || got.missing_value),
  if (iscellstr(f) || ischar(f)),
    ind = find(f == spval);
    f(ind) = spval;
  else
    ind = find(abs(f-spval) < 4*eps(double(spval)));
  end
end

%  Scale data and/or add offset value.

if (isnumeric(f)),
  if (got.add_offset),
    if (got.scale_factor),
      switch nctype
        case {nc_constant('nc_float'),                                  ...
              nc_constant('nc_double')}
          f = f.*scale+offset;
        case {nc_constant('nc_int'),                                    ...
              nc_constant('nc_short'),                                  ...
              nc_constant('nc_byte')}
          f = double(f).*scale+offset;
      end
    else
      switch nctype
        case {nc_constant('nc_float'),                                  ...
              nc_constant('nc_double')}
          f = f+offset;
        case {nc_constant('nc_int'),                                    ...
              nc_constant('nc_short'),                                  ...
              nc_constant('nc_byte')}
          f = double(f)+offset;
      end
    end
  elseif (got.scale_factor);
    switch nctype
      case {nc_constant('nc_float'),                                    ...
            nc_constant('nc_double')}
        f = f.*scale;
      case {nc_constant('nc_int'),                                      ...
            nc_constant('nc_short'),                                    ...
            nc_constant('nc_byte')}
        f = double(f).*scale;
    end
  end
end

%  Set fill values or missing values with specified ReplaceValue.

if (~isempty(ind) && isnumeric(f)),
  f(ind) = ReplaceValue;
end

%--------------------------------------------------------------------------
%  If desired, convert data to double precision.
%--------------------------------------------------------------------------

if (~PreserveType) && ~(iscellstr(f) || ischar(f)),
  f = double(f);
end

%--------------------------------------------------------------------------
% Close NetCDF file.
%--------------------------------------------------------------------------

status = mexnc('ncclose',ncid);
if (status == -1),
  error('NC_READ_MEXNC: ncclose - unable to close NetCDF file.')
end

return
