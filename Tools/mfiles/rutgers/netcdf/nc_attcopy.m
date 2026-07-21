function status = nc_attcopy(inpfile, outfile, varargin)

%
% NC_ATTCOPY:  copies a global or variable NetCDF attribute between files
%
% status = nc_attcopy(inpfile, outfile, Vname)
%
% This function copies a global or variable attribute from input to output
% NetCDF files. If the "Vname" argument is missing, it is assumed that
% "Aname" is a global attribute.
%
% On Input:
%
%    inpfile    Input  NetCDF filename (string)
%
%    outfile    Output NetCDF filename (string)
%
%    Vname      Variable name (string; optional)
%
% On Output:
%
%    status     Error flag
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Initialize.

Vname = [];
switch numel(varargin)
  case 1
    Vname = varargin{1};
end

% Initialize error flag.

status = 0;

if (isempty(Vname))
  got_var = false;
else
  got_var = true;
end

% Inquire about the contents of the input and output NetCDF files.

Iinp = nc_inq(inpfile);
Iout = nc_inq(outfile);

%--------------------------------------------------------------------------
% Add/modify a variable attribute.
%--------------------------------------------------------------------------

if (got_var)

% Get input NetCDF variable index and number of attributes.
 
  if (~any(strcmp({Iinp.Variables.Name}, Vname)))
    nc_inq(inpfile, true);
    disp(' ');
    error(['NC_ATTCOPY: cannot find input NetCDF variable: ', Vname]);
  else
    inpvar   = strcmp({Iinp.Variables.Name},  Vname);
    NinpAtts = length({Iinp.Variables(inpvar).Attributes.Name});
  end

% Get output NetCDF variable index and number of attributes.
  
  if (~any(strcmp({Iout.Variables.Name}, Vname)))
    nc_inq(outfile, true);
    disp(' ');
    error(['NC_ATTCOPY: cannot find output NetCDF variable: ', Vname]);
  else
    outvar   = strcmp({Iout.Variables.Name},  Vname);
    NoutAtts = length({Iout.Variables(outvar).Attributes.Name});
  end
  
% Copy variable attributes from input to output files.

  disp(' ');
  for i=1:NinpAtts
    Aname=Iinp.Variables(inpvar).Attributes(i).Name;
    Avalue=Iinp.Variables(inpvar).Attributes(i).Value;
    if (ischar(Avalue))
      disp(['Copying attribute: "',Aname,'" = ', Avalue]);
    else
      disp(['Copying attribute: "',Aname,'" = ', num2str(Avalue)]);
    end
    nc_attadd(outfile, Aname, Avalue, Vname);
  end

%--------------------------------------------------------------------------
%  Add/modify a global attribute.
%--------------------------------------------------------------------------
   
else

% Get input/output NetCDF global number of attributes.

  NinpAtts = length({Iinp.Attributes.Name});
  NoutAtts = length({Iout.Attributes.Name});

  disp(' ');
  for i=1:NinpAtts
    Aname=Iinp.Attributes(i).Name;
    Avalue=Iinp.Attributes(i).Value;
    if (ischar(Avalue))
      disp(['Copying global attribute: "',Aname,'" = ', Avalue]);
    else
      disp(['Copying global attribute: "',Aname,'" = ', num2str(Avalue)]);
    end
    nc_attadd(outfile, Aname, Avalue);
  end
  
end

return
