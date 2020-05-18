function [status]=nc_dfixed(ncfile_old, ncfile_new, Rname)
  
%
% NC_DFIXED:  Changes unlimited to fixed dimension in a NetCDF file
%
% status = nc_dfixed(ncfile_old, ncfile_new, Rname)
%
% This function creates a new NetCDF file with the requested record
% unlimited dimension changed to fixed dimension.
%
% On Input:
%
%    ncfile_old   Old NetCDF file name with unlimited dimension (string)
%
%    ncfile_new   New NetCDF file name with fixed dimension (string)
%
%    Rname        Record dimension name to change from unlimited to
%                   fixed (string)
%
% On Output:
%
%    status       Error flag
%

% svn $Id: nc_dfixed.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
  
%--------------------------------------------------------------------------
%  Get input file information structure.
%--------------------------------------------------------------------------

I = nc_inq(ncfile_old);

%--------------------------------------------------------------------------
%  Change requested dimension from unlimited to fixed (zero value).
%--------------------------------------------------------------------------

dnames = {I.Dimensions.Name};

if (~strcmp(dnames, Rname)),
  disp(' ');
  error(['NC_DFIXED: cannot find NetCDF dimension: ', Rname]);
else
  index = strcmp(dnames, Rname);
  if (~any(I.Dimensions(index).Unlimited)),
    disp(' ');
    error(['NC_DFIXED: requested dimension: ', Rname, ' is not unlimited']);
  end    
end

I.Dimensions(strcmp({I.Dimensions.Name}, Rname)).Unlimited=0;

%--------------------------------------------------------------------------
%  Inquire old file NetCDF format.
%--------------------------------------------------------------------------

ncid = netcdf.open(ncfile_old, 'NOWRITE');
frmt = netcdf.inqFormat(ncid);
netcdf.close (ncid);

%  Set new file creation mode.

mode = netcdf.getConstant('CLOBBER');               % overwite if existing
mode = bitor(mode,netcdf.getConstant(frmt));
 
%--------------------------------------------------------------------------
%  Create new NetCDF with fixed requested dimension.
%--------------------------------------------------------------------------

nc_create(ncfile_new, mode, I);

%--------------------------------------------------------------------------
%  Write out data into new NetCDF.
%--------------------------------------------------------------------------

vnames = {I.Variables.Name};

for value = vnames,

  field = char(value);
  F = nc_read(ncfile_old, field);

  ind=find(isnan(F));                           % overwrite NaNs with zeros
  if (~isempty(ind))
    F(ind)=0;
  end
  
  status = nc_write(ncfile_new, field, F);

end

return
