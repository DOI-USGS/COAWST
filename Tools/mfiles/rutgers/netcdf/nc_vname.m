function [vname,nvars]=nc_vname(fname)

%
% NC_VNAME:  Get the names of all variables in a NetCDF file
%
% [vname,nvars]=nc_vname(fname)
%
% This function gets the number and name of variables in requested
% NetCDF file.
%
% On Input:
%
%    fname      NetCDF file name or URL name (string)
%
% On Output:
%
%    vname      Variable names (padded string array)
%    nvars      Number of variables
%
% Warning:
% =======
%
% This function is OBSOLETE but it is kept for backward compatability.
% Use "nc_vnames" instead. The new funtion output arguments are different
% and the variable names are in a cell array and follow the native
% interface structure. It is a better way to manipulate string arrays.
% Notice that usage with the new "nc_vnames" is as follows:
%
% V=nc_vnames(fname)
% nvars=length(V.Variables); 
%
% for n=1:nvars,
%   name = char(V.Variables(n).Name);
%   switch name
%     ...
%   end
% end
%
  
% svn $Id: nc_vname.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Choose NetCDF file interface.

[method,~,~] = nc_interface(fname);

switch(method),
  case {'native'}
    [vname, nvars] = nc_vname_matlab(fname);   % Matlab native interface
  case {'java'}
    [vname, nvars] = nc_vname_java  (fname);   % SNCTOOLS JAVA interface
  case {'mexnc'}
    [vname, nvars] = nc_vname_mexnc (fname);   % MEXNC inteface
  otherwise
    error('NC_VNAME: unable to determine NetCDF processing interface');
end

return

%--------------------------------------------------------------------------

function [vname,nvars]=nc_vname_matlab(fname)

%
% NC_VNAME_MATLAB:  Get the names of all variables in a NetCDF file
%
% [vname,nvars]=nc_vname_matlab(fname)
%
% This function gets information about all the variable available in
% a regular or OpenDAP URL NetCDF file. It uses the native matlab
% function "ncinfo" that it is available in version 2012a or higher.
%

% Inquire information about requested variable.

Info = ncinfo(fname); 

% Extract variables information.

nvars = length(Info.Variables);

for n=1:nvars,
  name=Info.Variables(n).Name;
  lstr=length(name);
  vname(n,1:lstr)=name(1:lstr);
end

return

%--------------------------------------------------------------------------

function [vname,nvars]=nc_vname_java(fname)

%
% NC_VNAME_JAVA:  Get the names of all variables in a NetCDF file
%
% [vname,nvars]=nc_vname_java(fname)
%
% This function gets the number and name of variables in requested
% URL OpenDAP NetCDF file(s).
%

%  Inquire information from URL NetCDF file.

Info=nc_info(fname);

%  Extract requested variable information.

nvars=length(Info.Dataset);

for n=1:nvars,
  name=Info.Dataset(n).Name;
  lstr=length(name);
  vname(n,1:lstr)=name(1:lstr);
end

return

%--------------------------------------------------------------------------

function [vname,nvars]=nc_vname_mexnc(fname)

%
% NC_VNAME:  Get the names of all variables in a NetCDF file
%
% [vname,nvars]=nc_vname_mexnc(fname)
%
% This function gets the number and name of variables in requested
% NetCDF file.  It cannot process URL OpenDAP file(s).

%  Open NetCDF file.
 
[ncid]=mexnc('ncopen',fname,'nc_nowrite');
if (ncid == -1),
  error(['NC_VNAME: ncopen - unable to open file: ', fname]);
end
 
%  Inquire about the contents of NetCDf file. Display information.

[~,nvars,~,~,status]=mexnc('ncinquire',ncid);
if (status == -1),
  error([ 'NC_VNAME: ncinquire - error while inquiring file: ' fname]);
end,

%  Inquire information about all the variables.

n=0;
for i=0:nvars-1,
  [name,~,~,~,~,status]=mexnc('ncvarinq',ncid,i);
  if (status == -1),
    error(['NC_VNAME: ncvarinq - unable to inquire about variable ID: ',...
          num2str(i)]);
  else
    n=n+1;
    lstr=length(name);
    vname(n,1:lstr)=name(1:lstr);
  end
end

%  Close NetCDF file.

[status]=mexnc('ncclose',ncid);
if (status == -1),
  error('NC_VNAME: ncclose - unable to close NetCDF file.');
end

return
