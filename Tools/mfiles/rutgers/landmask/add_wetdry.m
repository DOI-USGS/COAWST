function add_wetdry(ncfile, mask_rho)

%
% ADD_WETDRY:  Adds wet/dry masks to a ROMS NetCDF file
%
% add_wetdry(ncfile, mask_rho)
%
% Adds wet/dry mask array to a ROMS NetCDF file.  If any of mask already
% exists, it will do nothing to define and process such variable.
%  
% This function can be use to attach the wet/dry mask at PSI-points. This
% fields was not defined in old versions of ROMS.  If wet/dry mask at
% RHO-points "mask_rho" is provided, it will recompute the other masks
% and will write those missing in the NetCDF.
%
% On Input:
%
%    ncfile      ROMS NetCDF file name (string)
%
%    mask_rho    Wet/dry mask at RHO-points (2D array; OPTIONAL)
%

% svn $Id: add_wetdry.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire grid NetCDF file about wet/dry mask variables.
%--------------------------------------------------------------------------

I=ncinfo(ncfile);
nvars=length(I.Variables);

% Set spherical switch.

index = strcmp({I.Variables.Name},'spherical');
if (any(index)),
  spherical=ncread(ncfile,'spherical');
  if (ischar(spherical)),
    if (spherical == 'T' || spherical == 't')
      spherical = true;
    else
      spherical = false;
    end
  end
else
  spherical=false;
end

% Set internal switches. 
  
if (nargin > 1),
  ComputeMask=true;
else
  ComputeMask=false;
end

got.mask_psi=false; define.mask_psi=true; Vname.mask_psi='wetdry_mask_psi';
got.mask_rho=false; define.mask_rho=true; Vname.mask_rho='wetdry_mask_rho';
got.mask_u  =false; define.mask_u  =true; Vname.mask_u  ='wetdry_mask_u';
got.mask_v  =false; define.mask_v  =true; Vname.mask_v  ='wetdry_mask_v';

for n=1:nvars,
  name=char(I.Variables(n).Name);
  switch name
    case {Vname.mask_psi}
      got.mask_psi=true;
      define.mask_psi=false;
    case {Vname.mask_rho}
      got.mask_rho=true;
      define.mask_rho=false;
    case {Vname.mask_u}
      got.mask_u=true;
      define.mask_u=false;
    case {Vname.mask_v}
      got.mask_v=true;
      define.mask_v=false;
  end
end

%--------------------------------------------------------------------------
%  If appropriate, define Land/Sea mask variables.
%--------------------------------------------------------------------------

if (define.mask_psi || define.mask_rho || define.mask_u || define.mask_v), 

% Open GRID NetCDF file and put into define mode.

  ncid=netcdf.open(ncfile,'WRITE');
  netcdf.reDef(ncid);

% Inquire about dimensions.

  Dname.xp='xi_psi';    Dname.yp='eta_psi';
  Dname.xr='xi_rho';    Dname.yr='eta_rho';
  Dname.xu='xi_u';      Dname.yu='eta_u';
  Dname.xv='xi_v';      Dname.yv='eta_v';

  ndims=length(I.Dimensions);

  for n=1:ndims,
    name=char(I.Dimensions(n).Name);
    switch name
      case {Dname.xp}
        Dsize.xp=I.Dimensions(n).Length;
        did.xp  =n-1;
      case {Dname.yp}
        Dsize.yp=I.Dimensions(n).Length;
        did.yp  =n-1;
      case {Dname.xr}
        Dsize.xr=I.Dimensions(n).Length;
        did.xr  =n-1;
      case {Dname.yr}
        Dsize.yr=I.Dimensions(n).Length;
        did.yr  =n-1;
      case {Dname.xu}
        Dsize.xu=I.Dimensions(n).Length;
        did.xu  =n-1;
      case {Dname.yu}
        Dsize.yu=I.Dimensions(n).Length;
        did.yu  =n-1;
      case {Dname.xv}
        Dsize.xv=I.Dimensions(n).Length;
        did.xv  =n-1;
      case {Dname.yv}
        Dsize.yv=I.Dimensions(n).Length;
        did.yv  =n-1;
    end
  end
  
%  Define wet/dry mask.

  if (define.mask_psi),
    varid=netcdf.defVar(ncid,Vname.mask_psi,                            ...
                        netcdf.getConstant('nc_double'),                ...
                        [did.xp did.yp]);
    netcdf.putAtt(ncid,varid,'long_name','wet/dry mask on PSI-points');
    netcdf.putAtt(ncid,varid,'flag_values',[0, 1]);
    netcdf.putAtt(ncid,varid,'flag_meanings','land water');
    if (spherical)
      netcdf.putAtt(ncid,varid,'coordinates','lon_psi lat_psi');
    else
      netcdf.putAtt(ncid,varid,'coordinates','x_psi y_psi');
    end
  end

  if (define.mask_rho),
    varid=netcdf.defVar(ncid,Vname.mask_rho,                            ...
                        netcdf.getConstant('nc_double'),                ...
                        [did.xr did.yr]);
    netcdf.putAtt(ncid,varid,'long_name','wet/dry mask on RHO-points');
    netcdf.putAtt(ncid,varid,'flag_values',[0, 1]);
    netcdf.putAtt(ncid,varid,'flag_meanings','land water');
    if (spherical)
      netcdf.putAtt(ncid,varid,'coordinates','lon_rho lat_rho');
    else
      netcdf.putAtt(ncid,varid,'coordinates','x_rho y_rho');
    end
  end

  if (define.mask_u),
    varid=netcdf.defVar(ncid,Vname.mask_u,                              ...
                        netcdf.getConstant('nc_double'),                ...
                        [did.xu did.yu]);
    netcdf.putAtt(ncid,varid,'long_name','wet/dry mask on U-points');
    netcdf.putAtt(ncid,varid,'flag_values',[0, 1]);
    netcdf.putAtt(ncid,varid,'flag_meanings','land water');
    if (spherical)
      netcdf.putAtt(ncid,varid,'coordinates','lon_u lat_u');
    else
      netcdf.putAtt(ncid,varid,'coordinates','x_u y_u');
    end
  end

  if (define.mask_v),
    varid=netcdf.defVar(ncid,Vname.mask_v,                              ...
                        netcdf.getConstant('nc_double'),                ...
                        [did.xv did.yv]);
    netcdf.putAtt(ncid,varid,'long_name','wet/dry mask on V-points');
    netcdf.putAtt(ncid,varid,'flag_values',[0, 1]);
    netcdf.putAtt(ncid,varid,'flag_meanings','land water');
    if (spherical)
      netcdf.putAtt(ncid,varid,'coordinates','lon_v lat_v');
    else
      netcdf.putAtt(ncid,varid,'coordinates','x_v y_v');
    end
  end

%  Leave definition mode and close NetCDF file.

  netcdf.endDef(ncid);
  netcdf.close(ncid);

end

%--------------------------------------------------------------------------
%  Compute and write out wet/dry mask(s) into NetCDF file.
%--------------------------------------------------------------------------

if (ComputeMask),
  [umask,vmask,pmask]=uvp_masks(mask_rho)

  if (define.mask_psi),
    ncwrite(ncfile,Vname.mask_psi,double(pmask));
  end

  if (define.mask_rho),
    ncwrite(ncfile,Vname.mask_rho,double(mask_rho));
  end

  if (define.mask_u),
    ncwrite(ncfile,Vname.mask_u,double(umask));
  end

  if (define.mask_v),
    ncwrite(ncfile,Vname.mask_psi,double(vmask));
  end
end

return
