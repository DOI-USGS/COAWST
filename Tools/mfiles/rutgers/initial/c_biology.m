function [status]=c_biology(S);

%
% C_BIOLOGY:  Defines biology model variables to a ROMS NetCDF file
%
% [status]=c_biology(S)
%
% This function defines biology model variables to an existing initial
% conditions or climatology ROMS NetCDF file.
%
% On Input:
%
%    S           Initial condidions creation parameters (structure array):
%
%                  S.ncname           NetCDF file name
%                  S.spherical        Spherical grid switch
%
% On Output:
%
%    status      Error flag.
%

% svn $Id: c_biology.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
  
%  Check required fields in structure S.

if (isfield(S,'ncname')),
  fname=S.ncname;
else,
  error([ 'C_BIOLOGY - cannot NetCDF file name field: ncname, ', ...
          'in structure array S']);
end,

%  Set variables type.

vartype=nc_constant('nc_double');

%  Inquire NetCDF file dimensions and their IDs.

D=nc_dinfo(fname);

for n=1:length(D),
  name=char(D(n).Name);
  switch name
    case 'xi_rho',
      did.xr=D(n).dimid;
      Lr=D(n).Length;
    case 'eta_rho',
      did.yr=D(n).dimid;
      Mr=D(n).Length;
    case 's_rho',
      did.Nr=D(n).dimid;
      Nr=D(n).Length;
    case 'ocean_time',
      did.time=D(n).dimid;
      rec=D(n).Length;
  end
end  

%  Horizontal grid variables.

Vname.time='ocean_time';

if (spherical),
  Vname.coord='lon_rho lat_rho s_rho ocean_time';
else,
  Vname.coord='x_rho y_rho s_rho ocean_time';
end,

%---------------------------------------------------------------------------
%  Open existing initial conditions or climatology file.
%---------------------------------------------------------------------------

%  Open NetCDF file.

[ncid]=mexnc('open',fname,'nc_write');
if (ncid < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BIOLOGY: open - unable to open file: ', fname]);
  return
end

%  Put open file into define mode.

[status]=mexnc('redef',ncid);
if (status < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  [status]=mexnc('close',ncid);
  error(['C_BIOLOGY: redef - unable to put in definition mode: ',fname]);
  return
end,

%---------------------------------------------------------------------------
%  Define variables according to biology model.
%---------------------------------------------------------------------------

switch S.biomodel,

  case {'BIO_FENNEL'}

    Var.name        = 'NO3';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'nitrate concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'NH4';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'ammonium concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'chlorophyll';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'chlorophyll concentration';
    Var.units       = 'milligrams meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
    
    Var.name        = 'phytoplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'phytoplankton biomass';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'zooplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'zooplankton biomass';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
    
    Var.name        = 'SdetritusN';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'small fraction nitrogen detritus concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'LdetritusN';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'large fraction nitrogen detritus concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'SdetritusC';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'small fraction carbon detritus concentration';
    Var.units       = 'millimole_carbon meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'LdetritusC';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'large fraction carbon detritus concentration';
    Var.units       = 'millimole_carbon meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
    
    Var.name        = 'TIC';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'total inorganic carbon';
    Var.units       = 'millimole_carbon meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'alkalinity';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'total alkalinity';
    Var.units       = 'milliequivalents meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'SdetritusN';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'small fraction nitrogen detritus concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'oxygen';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'dissolved oxygen concentration';
    Var.units       = 'millimole_oxygen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
    
 case {'NEMURO'}

    Var.name        = 'NO3';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'nitrate concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'NH4';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'ammonium concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'PON';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'particulate organic nitrogen concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'DON';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'dissolved organic nitrogen concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'SiOH4';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'silicate concentration';
    Var.units       = 'millimole_silica meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'opal';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'particulate organic silica concentration';
    Var.units       = 'millimole_silica meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'nanophytoplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'nanophytoplankton biomass';
    Var.units       = 'millimole meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
    
    Var.name        = 'diatom';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'diatom biomass';
    Var.units       = 'millimole meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'microzooplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'microzooplankton';
    Var.units       = 'millimole meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'mesozooplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'mesozooplankton';
    Var.units       = 'millimole meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'Pzooplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'predator-zooplankton biomass';
    Var.units       = 'millimole meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
    
  case {'NPZD_FRANKS', 'NPZD_POWELL'}
  
    Var.name        = 'NO3';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'nitrate concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'phytoplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'phytoplankton biomass';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'zooplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'zooplankton biomass';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'phytoplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'phytoplankton biomass';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'detritus';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'detritus concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

  case {'NPZD_IRON'}
  
    Var.name        = 'NO3';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'nitrate concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'iron';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'available dissolved iron concentration';
    Var.units       = 'millimole_iron meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
    
    Var.name        = 'phytoplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'phytoplankton biomass';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'phytoplanktonFe';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'phytoplankton, associated iron concentration';
    Var.units       = 'millimole_iron meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
    
    Var.name        = 'zooplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'zooplankton biomass';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'phytoplankton';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'phytoplankton biomass';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

    Var.name        = 'detritus';
    Var.type        = vartype;
    Var.dimid       = [did.time did.Nr did.yr did.xr];
    Var.long_name   = 'detritus concentration';
    Var.units       = 'millimole_nitrogen meter-3';
    Var.time        = Vname.time;
    Var.coordinates = Vname.coordinates;
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var

end,

%---------------------------------------------------------------------------
%  Close NetCDF file.
%---------------------------------------------------------------------------

%  Exit definition mode.

[status]=mexnc('enddef',ncid);
if (status < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BIOLOGY: enddef - unable to exit definition mode: ',fname]);
  return
end,

%  Close NetCDF file.

[cstatus]=mexnc('ncclose',ncid);
if (cstatus < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BIOLOGY: ncclose - unable to close NetCDF file: ', fname]);
end

return
