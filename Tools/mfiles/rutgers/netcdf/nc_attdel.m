function status = nc_attdel(ncfile, Aname, varargin)

%
% NC_ATTDEL:  Delete requested NetCDF attribute
%
% status = nc_attdel(ncfile, Aname, Vname);
%
% This function deletes requested global or variable attribute in a NetCDF
% file. If the "Vname" argument is missing, it is assumed that "Aname" is a
% global attribute.
%
% On Input:
%
%    ncfile      NetCDF file name (character string)
%
%    Aname      Attribute name (character string)
%
%    Vname      Variable name (character string; optional)
%
% On Output:
%
%    status     Error flag
%

% svn $Id: nc_attdel.m 711 2014-01-23 20:36:13Z arango $
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
    status = nc_attdel_matlab(ncfile, Aname, Vname, Info);
  case {'mexnc'}
    status = nc_attdel_mexnc (ncfile, Aname, Vname, Info);
  case {'java'}
    error('NC_ATTDEL: it is not possible to write into an OpenDAP file.');
  otherwise
    error('NC_ATTDEL: unable to determine NetCDF processing interface.');
end

return

%--------------------------------------------------------------------------

function status = nc_attdel_matlab(ncfile, Aname, Vname, Info)

%
% NC_ATTDEL_MATLAB:  Delete requested NetCDF attribute
%
% status = nc_attdel_matlab(ncfile, Aname, Vname, Info)
%
% This function deletes requested global or variable in a NetCDF file
% using the native Matlab interface. If the "Vname" argument is missing,
% it is assumed that "Aname" is a global attribute.
%
% On Input:
%
%    ncfile      NetCDF file name (string)
%
%    Aname      Attribute name (string)
%
%    Vname      Variable name (string)
%
%    Info       NetCDF information structure (struct array)
%
% On Output:
%
%    status     Error flag
%

%  Initialize error flag.

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
% Delete a variable attribute.
%---------------------------------------------------------------------------

if (got_var),

% Get variable ID.

  if (~any(strcmp({Info.Variables.Name}, Vname)))
    nc_inq(ncfile, true);
    disp(' ');
    error(['NC_ATTADD_MEXNC: cannot find NetCDF variable: ',Vname]);
  end
  varid = netcdf.inqVarID(ncid, Vname);
  ivar  = strcmp({Info.Variables.Name}, Vname);

% Delete requested variable attribute.

  if (any(strcmp({Info.Variables(ivar).Attributes.Name}, Aname))),
    netcdf.delAtt(ncid, varid, Aname);  
  else
    disp(' ');
    disp(['Requested attribute "',Aname,'" not found in variable "',    ...
          Vname,'".']);
  end
  
%---------------------------------------------------------------------------
% Delete a global attribute.
%---------------------------------------------------------------------------
   
else

% Delete requested global attribute.

  if (any(strcmp({Info.Attributes.Name}, Aname)))
    varid  = netcdf.getConstant('GLOBAL');
    netcdf.delAtt(ncid, varid, Aname);  
  else
    disp(' ');
    disp(['Requested global attribute "',Aname,'" not found in: ', ncfile]);
  end

end

% Exit definition mode and close NetCDF file.

netcdf.endDef(ncid);
netcdf.close(ncid);

return

%--------------------------------------------------------------------------

function status = nc_attdel_mexnc(ncfile, Aname, Vname, Info)

%
% NC_ATTDEL_MEXNC:  Delete requested NetCDF attribute
%
% status = nc_attdel_mexnc(ncfile, Aname, Vname, Info)
%
% This function deletes requested global or variable in a NetCDF file
% using the MEXNC interface. If the "Vname" argument is missing, it is
% assumed that "Aname" is a global attribute.
%
% On Input:
%
%    ncfile      NetCDF file name (string)
%
%    Aname      Attribute name (string)
%
%    Vname      Variable name (string)
%
%    Info       NetCDF information structure (struct array)
%
% On Output:
%
%    status     Error flag
%

%  Initialize error flag.

if (isempty(Vname)),
  got_var = false;
else
  got_var = true;
end

%  Open NetCDF file.

ncid = mexnc('open',ncfile,'nc_write');
if (ncid < 0),
  status = -1;
  disp(' ');
  error(['NC_ATTDEL_MEXNC: open - unable to open file: ', ncfile]);
end

%  Put open file into define mode.

status = mexnc('redef',ncid);
if (status < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  mexnc('close',ncid);
  error(['NC_ATTDEL_MEXNC: redef - unable to put in definition mode: ', ...
         ncfile]);
end

%---------------------------------------------------------------------------
%  Deleting a variable attribute.
%---------------------------------------------------------------------------

if (got_var),

%  Get variable ID.

  if (~any(strcmp({Info.Variables.Name}, Vname)))
    nc_inq(ncfile, true);
    disp(' ');
    error(['NC_ATTADD_MEXNC: cannot find NetCDF variable: ',Vname]);
  end
  varid = mexnc('inq_varid',ncid,Vname);

%  Inquire number of variable attributes.

  [nvatts,status] = mexnc('inq_varnatts',ncid,varid);
  if (status < 0),
    disp(' ');
    disp(mexnc('strerror',status));
    error(['NC_ATTDEL_MEXNC: inq_varnatts - unable to inquire number ', ...
	   'of variable attributes: ',Vname]);
  end

%  Delete requested variable attribute.

  found = false;
  
  for i=0:nvatts-1
   [attnam,status] = mexnc('inq_attname',ncid,varid,i);
   if (status < 0),
     disp(' ');
     disp(mexnc('strerror',status));
     error(['NC_ATTDEL_MEXNC: inq_attname: error while inquiring ',     ...
	    'attribute: ',num2str(i)]);
   end
   if (strcmp(Aname, attnam)),
     status = mexnc('del_att',ncid,varid,attnam);
     if (status < 0)
       disp(' ');
       disp(mexnc('strerror',status));
       error(['NC_ATTDEL_MEXNC: del_att - error while deleting ',       ...
              'attribute: ' Aname]);
     end,
     found = true;
     break
   end
 end
 if (~found),
   disp('  ');
   disp(['Requested attribute "',Aname,'" not found in variable "',     ...
         Vname,'".']);
 end
  
%---------------------------------------------------------------------------
% Deleting a global attribute.
%---------------------------------------------------------------------------
   
else

% Inquire number of global attributes.

  [natts,status] = mexnc('inq_natts',ncid);
  if (status < 0),
    disp(' ');
    disp(mexnc('strerror',status));
    error(['NC_ATTDEL_MEXNC: inq_natts - unable to inquire number ',    ...
	   'of global attributes: ',ncfile]);
  end
  
% Delete requested global attribute.

  found = false;
  
  for i=0:natts-1
    [attnam,status] = mexnc('inq_attname',ncid,nc_global,i);
    if (status < 0),
      disp(' ');
      disp(mexnc('strerror',status));
      error(['NC_ATTDEL_MEXNC: inq_attname: error while inquiring ',    ...
             'attribute: ', num2str(i)]);
    end
    if (strcmp(Aname, attnam)),
      status = mexnc('del_att',ncid,nc_global,attnam);
      if (status < 0)
        disp('  ');
        disp(mexnc('strerror',status));
        error(['NC_ATTDEL_MEXNC: del_att - error while deleting ',      ...
               'attribute: ' Aname]);
      end,
      found = true;
      break
    end
  end
  if (~found),
    disp(' ');
    disp(['Requested global attribute "',Aname,'" not found in: ', ncfile]);
  end
  
end

%  Exit definition mode.

status = mexnc('enddef',ncid);
if (status < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  error(['NC_ATTDEL_MEXNC: enddef - unable to exit definition mode: ',  ...
          ncfile]);
end

%  Close NetCDF file.

cstatus = mexnc('ncclose',ncid);
if (cstatus < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  error(['NC_ATTDEL_MEXNC: ncclose - unable to close NetCDF file: ',    ...
         ncfile]);
end

return
