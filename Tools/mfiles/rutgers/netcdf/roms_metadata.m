function V = roms_metadata(Vname,varargin)

%
% ROMS_METADATA:  Sets metadata structure for a ROMS NetCDF variable
%
% V = roms_metadata(Vname,spherical,nctype,Unlimited)
%
% This function sets the metadata structure, V, for requested ROMS
% NetCDF variable.  The structure contains the same fields that are
% returned by "nc_inq" or native interface "ncinfo".  This can be
% used elsewhere to create a NetCDF file in a compact way.
%
% On Input:
%
%    Vname      ROMS NetCDF variable to process (string)
%
%    spherical  Spherical switch (logical, OPTIONAL, default=true)
%
%    nctype     Variable datatype (string, OPTIONAL, default='nc_double')
%
%                 'nc_int'       integer variable
%                 'nc_float'     single precision variable
%                 'nc_double'    double precision variable
%
%    Unlimited  Variable has unlimited time dimension (logical, OPTIONAL,
%                                                      default=true)
%
% On Ouput:
%
%    V          Requested variable metadata structute (struct array):
%
%                 V.Name
%
%                 V.Dimensions(:).Name
%                 V.Dimensions(:).Length
%                 V.Dimensions(:).Unlimited
%
%                 V.Size
%
%                 V.Attributes(:).Name
%                 V.Attributes(:).Value
%
%                 V.Cgridtype.Name
%                 V.Cgridtype.Value
%
%                 V.Datatype
%
%                 V.ncType
%
% svn $Id: roms_metadata.m 722 2014-03-14 00:53:34Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMVariables.txt                   Hernan G. Arango      %
%=========================================================================%

% Set optional arguments.

spherical = true;
nctype    = 'nc_double'; Datatype  = 'double';
Unlimited = true;

switch numel(varargin)
 case 1
   spherical = varargin{1};
 case 2
   spherical = varargin{1};
   nctype    = varargin{2};
 case 3
   spherical = varargin{1};
   nctype    = varargin{2};
   Unlimited = varargin{3};
end

if (numel(varargin) > 1),
 if ~(isempty(strfind(nctype,'_')))
   is=strfind(nctype,'_')+1;
   ie=length(nctype);
   Datatype=lower(nctype(is:ie));
 else
   Datatype=lower(nctype);
 end
end

%==========================================================================
%  Set NetCDF Rmetadata structure for requested ROMS variable.
%==========================================================================

switch Vname
 
%--------------------------------------------------------------------------
%  Grid variables.
%--------------------------------------------------------------------------
 
 case 'spherical'
    V.Name                    = Vname;
    V.Dimensions              = [];
    V.Size                    = 1;
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'grid type logical switch';
    V.Attributes(2).Name      = 'flag_values';
    V.Attributes(2).Value     = [int32(0) int32(1)];
    V.Attributes(3).Name      = 'flag_meanings';
    V.Attributes(3).Value     = 'Cartesian spherical';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'int32';
    V.ncType                  = nc_constant('nc_int');
  case 'Vtransform'
    V.Name                    = Vname;
    V.Dimensions              = [];
    V.Size                    = 1;
    V.Attributes.Name         = 'long_name';
    V.Attributes.Value        = 'vertical terrain-following transformation equation';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'int32';
    V.ncType                  = nc_constant('nc_int');
  case 'Vstretching'
    V.Name                    = Vname;
    V.Dimensions              = [];
    V.Size                    = 1;
    V.Attributes.Name         = 'long_name';
    V.Attributes.Value        = 'vertical terrain-following stretching function';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'int32';
    V.ncType                  = nc_constant('nc_int');
  case 'theta_s'
    V.Name                    = Vname;
    V.Dimensions              = [];
    V.Size                    = 1;
    V.Attributes.Name         = 'long_name';
    V.Attributes.Value        = 'S-coordinate surface control parameter';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'theta_b'
    V.Name                    = Vname;
    V.Dimensions              = [];
    V.Size                    = 1;
    V.Attributes.Name         = 'long_name';
    V.Attributes.Value        = 'S-coordinate bottom control parameter';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'Tcline'
    V.Name                    = Vname;
    V.Dimensions              = [];
    V.Size                    = 1;
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'S-coordinate surface/bottom layer width';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'hc'
    V.Name                    = Vname;
    V.Dimensions              = [];
    V.Size                    = 1;
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'S-coordinate parameter, critical depth';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 's_rho'
    V.Name                    = Vname;
    V.Dimensions.Name         = 's_rho';
    V.Dimensions.Length       = [];
    V.Dimensions.Unlimited    = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'S-coordinate at RHO-points';
    V.Attributes(2).Name      = 'valid_min';
    V.Attributes(2).Value     = -1;
    V.Attributes(3).Name      = 'valid_max';
    V.Attributes(3).Value     = 0;
    V.Attributes(4).Name      = 'positive';
    V.Attributes(4).Value     = 'up';
    V.Attributes(5).Name      = 'standard_name';
    V.Attributes(5).Value     = 'ocean_s_coordinate_g2';
    V.Attributes(6).Name      = 'formula_terms';
    V.Attributes(6).Value     = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 's_w'
    V.Name                    = Vname;
    V.Dimensions.Name         = 's_w';
    V.Dimensions.Length       = [];
    V.Dimensions.Unlimited    = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'S-coordinate at W-points';
    V.Attributes(2).Name      = 'valid_min';
    V.Attributes(2).Value     = -1;
    V.Attributes(3).Name      = 'valid_max';
    V.Attributes(3).Value     = 0;
    V.Attributes(4).Name      = 'positive';
    V.Attributes(4).Value     = 'up';
    V.Attributes(5).Name      = 'standard_name';
    V.Attributes(5).Value     = 'ocean_s_coordinate_g2';
    V.Attributes(6).Name      = 'formula_terms';
    V.Attributes(6).Value     = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'Cs_r'
    V.Name                    = Vname;
    V.Dimensions.Name         = 's_rho';
    V.Dimensions.Length       = [];
    V.Dimensions.Unlimited    = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'S-coordinate stretching curves at RHO-points';
    V.Attributes(2).Name      = 'valid_min';
    V.Attributes(2).Value     = -1;
    V.Attributes(3).Name      = 'valid_max';
    V.Attributes(3).Value     = 0;
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'Cs_w'
    V.Name                    = Vname;
    V.Dimensions.Name         = 's_w';
    V.Dimensions.Length       = [];
    V.Dimensions.Unlimited    = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'S-coordinate stretching curves at W-points';
    V.Attributes(2).Name      = 'valid_min';
    V.Attributes(2).Value     = -1;
    V.Attributes(3).Name      = 'valid_max';
    V.Attributes(3).Value     = 0;
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'h'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'bathymetry at RHO-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho';
    else
      V.Attributes(3).Value   = 'x_rho y_rho';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'f'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'Coriolis parameter at RHO-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'second-1';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho';
    else
      V.Attributes(3).Value   = 'x_rho y_rho';
    end      
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'pm'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'curvilinear coordinate metric in XI';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter-1';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho';
    else
      V.Attributes(3).Value   = 'x_rho y_rho';
    end      
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'pn'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'curvilinear coordinate metric in ETA';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter-1';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho';
    else
      V.Attributes(3).Value   = 'x_rho y_rho';
    end      
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'x_rho'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'x-locations of RHO-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'y_rho'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'y-locations of RHO-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'x_psi'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_psi';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_psi';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'x-locations of PSI-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Cgridtype.Name          = 'streamfunction point';
    V.Cgridtype.Value         = 2;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'y_psi'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_psi';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_psi';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'y-locations of PSI-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Cgridtype.Name          = 'streamfunction point';
    V.Cgridtype.Value         = 2;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'x_u'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_u';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'x-locations of U-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'y_u'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_u';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'y-locations of U-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'x_v'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_v';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'x-locations of V-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'y_v'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_v';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'y-locations of V-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'longitude';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'latitude';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_rho'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'longitude of RHO-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_rho'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'latitude of RHO-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_psi'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_psi';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_psi';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'longitude of PSI-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'streamfunction point';
    V.Cgridtype.Value         = 2;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_psi'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_psi';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_psi';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'latitude of PSI-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'streamfunction point';
    V.Cgridtype.Value         = 2;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_u'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_u';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'longitude of U-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_u'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_u';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'latitude of U-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_v'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_v';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'longitude of V-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_v'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_v';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'latitude of V-points';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'mask_rho'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'mask on RHO-points';
    V.Attributes(2).Name      = 'flag_values';
    V.Attributes(2).Value     = [double(0) double(1)];
    V.Attributes(3).Name      = 'flag_meanings';
    V.Attributes(3).Value     = 'land water';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho';
    else
      V.Attributes(4).Value   = 'x_rho y_rho';
    end      
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'mask_psi'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_psi';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_psi';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'mask on PSI-points';
    V.Attributes(2).Name      = 'flag_values';
    V.Attributes(2).Value     = [double(0) double(1)];
    V.Attributes(3).Name      = 'flag_meanings';
    V.Attributes(3).Value     = 'land water';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_psi lat_psi';
    else
      V.Attributes(4).Value   = 'x_psi y_psi';
    end      
    V.Cgridtype.Name          = 'streamfunction point';
    V.Cgridtype.Value         = 2;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'mask_u'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_u';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'mask on U-points';
    V.Attributes(2).Name      = 'flag_values';
    V.Attributes(2).Value     = [double(0) double(1)];
    V.Attributes(3).Name      = 'flag_meanings';
    V.Attributes(3).Value     = 'land water';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u lat_u';
    else
      V.Attributes(4).Value   = 'x_u y_u';
    end      
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'mask_v'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_v';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'mask on V-points';
    V.Attributes(2).Name      = 'flag_values';
    V.Attributes(2).Value     = [double(0) double(1)];
    V.Attributes(3).Name      = 'flag_meanings';
    V.Attributes(3).Value     = 'land water';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_v lat_v';
    else
      V.Attributes(4).Value   = 'x_v y_v';
    end      
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'visc_factor'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'horizontal viscosity sponge factor';
    V.Attributes(2).Name      = 'valid_min';
    V.Attributes(2).Value     = double(0);
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho';
    else
      V.Attributes(3).Value   = 'x_rho y_rho';
    end      
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'diff_factor'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'horizontal diffusivity sponge factor';
    V.Attributes(2).Name      = 'valid_min';
    V.Attributes(2).Value     = double(0);
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho';
    else
      V.Attributes(3).Value   = 'x_rho y_rho';
    end      
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
    
%--------------------------------------------------------------------------
%  Boundary conditions grid variables.
%--------------------------------------------------------------------------

  case 'lon_rho_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of RHO-points, ',           ...
                                 'western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_rho_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of RHO-points, ',            ...
                                 'western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_rho_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of RHO-points, ',           ...
                                 'eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_rho_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of RHO-points, ',            ...
                                 'eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_rho_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of RHO-points, ',           ...
                                 'southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_rho_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of RHO-points, ',            ...
                                 'southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_rho_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of RHO-points, ',           ...
                                 'northern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_rho_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of RHO-points, ',            ...
                                 'northern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_u_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of U-points, ',             ...
                                 'western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_u_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of U-points, ',              ...
                                 'western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_u_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of U-points, ',             ...
                                 'eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_u_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of U-points, ',              ...
                                 'eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_u_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of U-points, ',             ...
                                 'southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_u_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of U-points, ',              ...
                                 'southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_u_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of U-points, ',             ...
                                 'northern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_u_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of U-points, ',              ...
                                 'northern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_v_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of V-points, ',             ...
                                 'western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_v_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of V-points, ',              ...
                                 'western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_v_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of V-points, ',             ...
                                 'eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_v_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of V-points, ',              ...
                                 'eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_v_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of V-points, ',             ...
                                 'southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_v_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of V-points, ',              ...
                                 'southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lon_v_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['longitude of V-points, ',             ...
                                 'northern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_east';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'longitude';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lat_v_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['latitude of V-points, ',              ...
                                 'northern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degree_north';
    V.Attributes(3).Name      = 'standard_name';
    V.Attributes(3).Value     = 'latitude';
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');

%--------------------------------------------------------------------------
%  Surface atmospheric forcing variables.
%--------------------------------------------------------------------------

  case 'cloud_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'cloud_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'cloud fraction time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'cloud'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'cloud_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'cloud fraction';
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'cloud_time';
    V.Attributes(3).Name      = 'coordinates';
    V.Attributes(3).Value     = 'lon lat cloud_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'pair_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'pair_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface air pressure time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'Pair'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'pair_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface air pressure';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millibar';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'pair_time';
    V.Attributes(4).Name      = 'coordinates';
    V.Attributes(4).Value     = 'lon lat pair_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'qair_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'qair_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface air humidity time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'Qair'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'qair_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface air relative humidity';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'percentage';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'qair_time';
    V.Attributes(4).Name      = 'coordinates';
    V.Attributes(4).Value     = 'lon lat qair_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'tair_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'tair_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface air temperature time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'Tair'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'tair_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface air temperature';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Celsius';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'tair_time';
    V.Attributes(4).Name      = 'coordinates';
    V.Attributes(4).Value     = 'lon lat tair_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'lhf_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lhf_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'latent heat flux time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'latent'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'lhf_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'net latent heat flux';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Watt meter-2';
    V.Attributes(3).Name      = 'positive_value';
    V.Attributes(3).Value     = 'downward flux, heating';
    V.Attributes(4).Name      = 'negative_value';
    V.Attributes(4).Value     = 'upward flux, cooling';
    V.Attributes(5).Name      = 'time';
    V.Attributes(5).Value     = 'lrf_time';
    V.Attributes(6).Name      = 'coordinates';
    V.Attributes(6).Value     = 'lon lat lhf_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'lrf_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lrf_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'longwave radiation flux time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'lwrad'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'lrf_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'net longwave radiation flux';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Watt meter-2';
    V.Attributes(3).Name      = 'positive_value';
    V.Attributes(3).Value     = 'downward flux, heating';
    V.Attributes(4).Name      = 'negative_value';
    V.Attributes(4).Value     = 'upward flux, cooling';
    V.Attributes(5).Name      = 'time';
    V.Attributes(5).Value     = 'lrf_time';
    V.Attributes(6).Name      = 'coordinates';
    V.Attributes(6).Value     = 'lon lat lrf_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'lwrad_down'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'lrf_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'downwelling longwave radiation flux';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Watt meter-2';
    V.Attributes(3).Name      = 'positive_value';
    V.Attributes(3).Value     = 'downward flux, heating';
    V.Attributes(4).Name      = 'negative_value';
    V.Attributes(4).Value     = 'upward flux, cooling';
    V.Attributes(5).Name      = 'time';
    V.Attributes(5).Value     = 'lrf_time';
    V.Attributes(6).Name      = 'coordinates';
    V.Attributes(6).Value     = 'lon lat lrf_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'sen_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'sen_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'sensible heat flux time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'sensible'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'sen_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'net sensible heat flux';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Watt meter-2';
    V.Attributes(3).Name      = 'positive_value';
    V.Attributes(3).Value     = 'downward flux, heating';
    V.Attributes(4).Name      = 'negative_value';
    V.Attributes(4).Value     = 'upward flux, cooling';
    V.Attributes(5).Name      = 'time';
    V.Attributes(5).Value     = 'lrf_time';
    V.Attributes(6).Name      = 'coordinates';
    V.Attributes(6).Value     = 'lon lat sen_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'shf_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'shf_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface net heat flux time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'shflux'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'shf_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface net heat flux';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Watt meter-2';
    V.Attributes(3).Name      = 'positive_value';
    V.Attributes(3).Value     = 'downward flux, heating';
    V.Attributes(4).Name      = 'negative_value';
    V.Attributes(4).Value     = 'upward flux, cooling';
    V.Attributes(5).Name      = 'time';
    V.Attributes(5).Value     = 'shf_time';
    V.Attributes(6).Name      = 'coordinates';
    if (spherical),
      V.Attributes(6).Value   = 'lon_rho lat_rho shf_time';
    else
      V.Attributes(6).Value   = 'x_rho y_rho shf_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'sms_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'sms_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface momentum stress time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'sustr'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_u';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'sms_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface u-momentum stress';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Newton meter-2';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'sms_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u lat_u sms_time';
    else
      V.Attributes(4).Value   = 'x_u y_u sms_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'svstr'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_v';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'sms_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface v-momentum stress';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Newton meter-2';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'sms_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_v lat_v sms_time';
    else
      V.Attributes(4).Value   = 'x_v y_v sms_time';
    end
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'swf_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'swf_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface net freshwater flux time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'swflux'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'swf_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface net freshwater flux (E-P)';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'centimeter day-1';
    V.Attributes(3).Name      = 'positive_value';
    V.Attributes(3).Value     = 'net evaporation';
    V.Attributes(4).Name      = 'negative_value';
    V.Attributes(4).Value     = 'net precipitation';
    V.Attributes(5).Name      = 'time';
    V.Attributes(5).Value     = 'swf_time';
    V.Attributes(6).Name      = 'coordinates';
    if (spherical),
      V.Attributes(6).Value   = 'lon_rho lat_rho swf_time';
    else
      V.Attributes(6).Value   = 'x_rho y_rho swf_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'srf_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'srf_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'shortwave radiation flux time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'swrad'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'srf_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'solar shortwave radiation flux';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Watt meter-2';
    V.Attributes(3).Name      = 'positive_value';
    V.Attributes(3).Value     = 'downward flux, heating';
    V.Attributes(4).Name      = 'negative_value';
    V.Attributes(4).Value     = 'upward flux, cooling';
    V.Attributes(5).Name      = 'time';
    V.Attributes(5).Value     = 'srf_time';
    V.Attributes(6).Name      = 'coordinates';
    V.Attributes(6).Value     = 'lon lat srf_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'rain_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'rain_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'rain fall time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'rain'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'rain_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'rain fall';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'kilogram meter-2 second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'rain_time';
    V.Attributes(4).Name      = 'coordinates';
    V.Attributes(4).Value     = 'lon lat rain_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'wind_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'wind_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface wind time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'Uwind'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wind_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface u-wind component';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'wind_time';
    V.Attributes(4).Name      = 'coordinates';
    V.Attributes(4).Value     = 'lon lat wind_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
 case 'Vwind'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wind_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface v-wind component';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'wind_time';
    V.Attributes(4).Name      = 'coordinates';
    V.Attributes(4).Value     = 'lon lat wind_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'sst_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'sst_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'sea surface temperature time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'SST'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'sst_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'sea surface temperature';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Celsius';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'sst_time';
    V.Attributes(4).Name      = 'coordinates';
    V.Attributes(4).Value     = 'lon lat sst_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'dQdSST'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'lon';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'lat';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'sst_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'surface net heat flux sensitivity to SST';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Watt meter-2 Celsius-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'sst_time';
    V.Attributes(4).Name      = 'coordinates';
    V.Attributes(4).Value     = 'lon lat sst_time';
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
    
%--------------------------------------------------------------------------
%  Wind-induced waves forcing variables.
%--------------------------------------------------------------------------

  case 'Dwave'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wave_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'wind-induced wave direction';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degrees';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'wave_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho wave_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho wave_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'Hwave'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wave_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'wind-induced significant wave height';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'wave_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho wave_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho wave_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'Lwave'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wave_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'wind-induced mean wavelength';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'wave_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho wave_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho wave_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'Pwave_top'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wave_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'wind-induced surface wave period';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'second';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'wave_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho wave_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho wave_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'Pwave_bot'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wave_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'wind-induced bottom wave period';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'second';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'wave_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho wave_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho wave_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'Ub_swan'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wave_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'SWAN, wind-induced bottom orbital velocity';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'wave_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho wave_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho wave_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'Ub_wave'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wave_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'wave model, wind-induced bottom orbital velocity';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'wave_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho wave_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho wave_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'Wave_dissip'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wave_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'wave dissipation';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Watts meter-2';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'wave_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho wave_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho wave_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'Wave_break'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'wave_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'percent wave breaking';
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'wave_time';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho wave_time';
    else
      V.Attributes(3).Value   = 'x_rho y_rho wave_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'wave_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'wave_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'wind-induce wave time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
    
%--------------------------------------------------------------------------
%  Tidal forcing variables.
%--------------------------------------------------------------------------

  case 'tide_period'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'tide_period';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'tide angular period';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'hour';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'tide_Ephase'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'tide_period';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'tidal elevation phase angle';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degrees, time of maximum elevation with respect chosen time origin';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho tide_period';
    else
      V.Attributes(3).Value   = 'x_rho y_rho tide_period';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'tide_Eamp'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'tide_period';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'tidal elevation amplitude';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho tide_period';
    else
      V.Attributes(3).Value   = 'x_rho y_rho tide_period';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'tide_Cphase'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'tide_period';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'tidal current phase angle';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'degrees, time of maximum velocity with respect chosen time origin';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho tide_period';
    else
      V.Attributes(3).Value   = 'x_rho y_rho tide_period';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'tide_Cangle'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'tide_period';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'tidal current inclination angle';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho tide_period';
    else
      V.Attributes(3).Value   = 'x_rho y_rho tide_period';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'tide_Cmin'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'tide_period';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'minimum tidal current, ellipse semi-minor axis';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho tide_period';
    else
      V.Attributes(3).Value   = 'x_rho y_rho tide_period';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'tide_Cmax'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'tide_period';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'maximum tidal current, ellipse semi-major axis';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho tide_period';
    else
      V.Attributes(3).Value   = 'x_rho y_rho tide_period';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

%--------------------------------------------------------------------------
%  River runoff.
%--------------------------------------------------------------------------
 
  case 'river_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'river_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'river runoff time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'river'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'river';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'river runoff identification number';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'river_Xposition'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'river';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'river XI-position at RHO-points';
    V.Attributes(2).Name      = 'valid_min';
    V.Attributes(2).Value     = 1;
    V.Attributes(3).Name      = 'valid_max';
    V.Attributes(3).Value     = [];
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'river_Eposition'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'river';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'river ETA-position at RHO-points';
    V.Attributes(2).Name      = 'valid_min';
    V.Attributes(2).Value     = 1;
    V.Attributes(3).Name      = 'valid_max';
    V.Attributes(3).Value     = [];
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'river_direction'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'river';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'river runoff direction';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'river_Vshape'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 's_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'river';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'river runoff mass transport vertical profile';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'river_transport'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'river';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'river_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'river runoff vertically integrated mass transport';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter3 second-1';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'river_temp'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'river';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(1).Name      = 's_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(3).Name      = 'river_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'river runoff potential temperature';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Celsius';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'river_time';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'river_salt'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'river';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(1).Name      = 's_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(3).Name      = 'river_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'river runoff salinity';
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'river_time';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');

%--------------------------------------------------------------------------
%  Physical state variables.
%--------------------------------------------------------------------------

  case 'clm_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'clm_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'climatology time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'ocean_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'ocean_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'time since intialization';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'second';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');
  case 'AKs'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_w';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'salinity vertical diffusion coefficient';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter2 second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_w ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_w ocean_time';
    end
    V.Cgridtype.Name          = 'w-velocity point';
    V.Cgridtype.Value         = 5;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'AKt'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_w';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'temperature vertical diffusion coefficient';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter2 second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_w ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_w ocean_time';
    end
    V.Cgridtype.Name          = 'w-velocity point';
    V.Cgridtype.Value         = 5;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'AKv'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_w';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'vertical viscosity coefficient';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter2 second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_w ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_w ocean_time';
    end
    V.Cgridtype.Name          = 'w-velocity point';
    V.Cgridtype.Value         = 5;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'rho'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'density anomaly';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'kilogram meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'omega'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_w';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'S-coordinate vertical momentum component';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_w ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_w ocean_time';
    end
    V.Cgridtype.Name          = 'w-velocity point';
    V.Cgridtype.Value         = 5;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'temp'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'potential temperature';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Celsius';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'salt'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'salinity';
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'ocean_time';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(3).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'u'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_u';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'u-momentum component';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u lat_u s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_u y_u s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'ubar'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_u';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'ocean_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'vertically integrated u-momentum component';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u lat_u ocean_time';
    else
      V.Attributes(4).Value   = 'x_u y_u ocean_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'v'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_v';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'v-momentum component';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_v lat_v s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_v y_v s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'vbar'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_v';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'ocean_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'vertically integrated v-momentum component';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_v lat_v ocean_time';
    else
      V.Attributes(4).Value   = 'x_v y_v ocean_time';
    end
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'w'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_w';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'vertical momentum component';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_w ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_w ocean_time';
    end
    V.Cgridtype.Name          = 'w-velocity point';
    V.Cgridtype.Value         = 5;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'zeta'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'ocean_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'free-surface';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

%--------------------------------------------------------------------------
%  Physical state boundary conditions variables.
%--------------------------------------------------------------------------

  case 'bry_time'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'bry_time';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'boundary conditions time';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'second';
    V.Cgridtype.Name          = 'none';
    V.Cgridtype.Value         = 0;
    V.Datatype                = 'double';
    V.ncType                  = nc_constant('nc_double');

  case 'temp_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['potential temperature, ',             ...
                                 'western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Celsius';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho_west lat_rho_west s_rho bry_time';
    else
      V.Attributes(4).Value   = 'x_rho_west y_rho_west s_rho bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'temp_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['potential temperature, ',             ...
                                 'eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Celsius';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho_east lat_rho_east s_rho bry_time';
    else
      V.Attributes(4).Value   = 'x_rho_east y_rho_east s_rho bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'temp_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['potential temperature, ',             ...
                                 'southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Celsius';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho_south lat_rho_south s_rho bry_time';
    else
      V.Attributes(4).Value   = 'x_rho_south y_rho_south s_rho bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'temp_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['potential temperature, ',             ...
                                 'northern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Celsius';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho_north lat_rho_north s_rho bry_time';
    else
      V.Attributes(4).Value   = 'x_rho_north y_rho_north s_rho bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'salt_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'salinity, western boundary condition';
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'bry_time';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho_west lat_rho_west s_rho bry_time';
    else
      V.Attributes(3).Value   = 'x_rho_west y_rho_west s_rho bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'salt_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'salinity, eastern boundary condition';
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'bry_time';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho_east lat_rho_east s_rho bry_time';
    else
      V.Attributes(3).Value   = 'x_rho_east y_rho_east s_rho bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'salt_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'salinity, southern boundary condition';
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'bry_time';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho_south lat_rho_south s_rho bry_time';
    else
      V.Attributes(3).Value   = 'x_rho_south y_rho_south s_rho bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'salt_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'salinity, northern boundary condition';
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'bry_time';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho_north lat_rho_north s_rho bry_time';
    else
      V.Attributes(3).Value   = 'x_rho_north y_rho_north s_rho bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'u_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['u-momentum component, ',              ...
                                 'western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u_west lat_u_west s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_u_west y_u_west s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'u_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['u-momentum component, ',              ...
                                 'eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u_east lat_u_east s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_u_east y_u_east s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'u_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['u-momentum component, ',              ...
                                 'southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u_south lat_u_south s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_u_south y_u_south s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'u_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['u-momentum component, ',              ...
                                 'northern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u_north lat_u_north s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_u_north y_u_north s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'ubar_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['vertically integrated u-momentum ',   ...
                                 'component, western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u_west lat_u_west bry_time';
    else
      V.Attributes(4).Value   = 'x_u_west y_u_west bry_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'ubar_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['vertically integrated u-momentum ',   ...
                                 'component, eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u_east lat_u_east bry_time';
    else
      V.Attributes(4).Value   = 'x_u_east y_u_east bry_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'ubar_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['vertically integrated u-momentum ',   ...
                                 'component, southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u_south lat_u_south bry_time';
    else
      V.Attributes(4).Value   = 'x_u_south y_u_south bry_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'ubar_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_u';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['vertically integrated u-momentum ',   ...
                                 'component, northern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_u_north lat_u_north bry_time';
    else
      V.Attributes(4).Value   = 'x_u_north y_u_north bry_time';
    end
    V.Cgridtype.Name          = 'u-velocity point';
    V.Cgridtype.Value         = 3;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'v_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['v-momentum component, ',              ...
                                 'western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_v_west lat_v_west s_rho bry_time';
    else
      V.Attributes(4).Value   = 'x_v_west y_v_west s_rho bry_time';
    end
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'v_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['v-momentum component, ',              ...
                                 'eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_v_east lat_v_east s_rho bry_time';
    else
      V.Attributes(4).Value   = 'x_v_east y_v_east s_rho bry_time';
    end
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'v_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 's_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'bry_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['v-momentum component, ',              ...
                                 'southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_v_south lat_v_south s_rho bry_time';
    else
      V.Attributes(4).Value   = 'x_v_south y_v_south s_rho bry_time';
    end
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'vbar_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['vertically integrated v-momentum ',   ...
                                 'component, western boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_v_west lat_v_west bry_time';
    else
      V.Attributes(4).Value   = 'x_v_west y_v_west bry_time';
    end
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'vbar_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['vertically integrated v-momentum ',   ...
                                 'component, eastern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_v_east lat_v_east bry_time';
    else
      V.Attributes(4).Value   = 'x_v_east y_v_east bry_time';
    end
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'vbar_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_v';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['vertically integrated v-momentum ',   ...
                                 'component, southern boundary condition'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_v_south lat_v_south bry_time';
    else
      V.Attributes(4).Value   = 'x_v_south y_v_south bry_time';
    end
    V.Cgridtype.Name          = 'v-velocity point';
    V.Cgridtype.Value         = 4;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'zeta_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'free-surface, eastern boundary condition';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho_west lat_rho_west bry_time';
    else
      V.Attributes(4).Value   = 'x_rho_west y_rho_west ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'zeta_west'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'free-surface, western boundary condition';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho_west lat_rho_west bry_time';
    else
      V.Attributes(4).Value   = 'x_rho_west y_rho_west bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'zeta_east'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'eta_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'free-surface, eastern boundary condition';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho_east lat_rho_east bry_time';
    else
      V.Attributes(4).Value   = 'x_rho_east y_rho_east bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'zeta_south'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'free-surface, southern boundary condition';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho_south lat_rho_south bry_time';
    else
      V.Attributes(4).Value   = 'x_rho_south y_rho_south bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'zeta_north'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'bry_time';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'free-surface, northern boundary condition';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'bry_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho_north lat_rho_north bry_time';
    else
      V.Attributes(4).Value   = 'x_rho_north y_rho_north bry_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
    
%--------------------------------------------------------------------------
%  Biology state variables.
%--------------------------------------------------------------------------

  case 'alkalinity'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'total alkalinity';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'milliequivalents meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'chlorophyll'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'chlorophyll concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_chlorophyll meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'detritus'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'detritus concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'diatom'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'diatom biomass';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'DON'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'dissolved organic nitrogen concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'iron'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'available dissolved iron concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'LdetritusC'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'large fraction carbon detritus concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_carbon meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'LdetritusN'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'large fraction nitrogen detritus concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'mesozooplankton'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'mesozooplankton biomass';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'microzooplankton'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'microzooplankton biomass';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'nanophytoplankton'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'nanophytoplankton biomass';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'NH4'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'ammonium concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'NO3'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'nitrate concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'opal'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'particulate organic silica concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_silica meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'oxygen'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'dissolved oxygen concentration functional';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_oxygen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'phytoplankton'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'phytoplankton concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'phytoplanktonFe'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'phytoplankton, associated iron concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_iron meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'PON'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'particulate organic nitrogen concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'Pzooplankton'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'predator-zooplankton biomass';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'SdetritusC'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'small fraction carbon detritus concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_carbon meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'SdetritusN'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'small fraction nitrogen detritus concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'SiOH4'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'silicate concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_silica meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype); 
  case 'zooplankton'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'zooplankton concentration';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_nitrogen meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
  case 'TIC'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'total inorganic carbon';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'millimole_carbon meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

%--------------------------------------------------------------------------
%  Sediment state variables.
%--------------------------------------------------------------------------

  case 'bed_age'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'Nbed';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'sediment bed layer age';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'second';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho Nbed ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho Nbed ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case 'bed_biodiff'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'Nbed';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'biodiffusivity at bottom of each layer';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter2 second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho Nbed ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho Nbed ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case 'bed_porosity'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'Nbed';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'sediment bed layer porosity';
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'ocean_time';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho Nbed ocean_time';
    else
      V.Attributes(3).Value   = 'x_rho y_rho Nbed ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case 'bed_thickness'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'Nbed';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'sediment bed layer thickness';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho Nbed ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho Nbed ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
 
  case 'erosion_stress'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'ocean_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'sediment median critical erosion stress';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'Newton meter-2';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case 'grain_density'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'ocean_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'sediment median grain density';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'kilogram meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case 'grain_diameter'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'ocean_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'sediment median grain diameter size';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
   
  case {'mud_01', 'mud_02', 'mud_03', 'mud_04',                         ...
        'mud_05', 'mud_06', 'mud_07', 'mud_08'}
    class = textscan(Vname, 'mud_ %d');
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['suspended cohesive sediment, ',       ...
                                 'size class ',                         ...
                                 num2str(class{:}, '%2.2i')];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'kilogram meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case {'mudfrac_01', 'mudfrac_02', 'mudfrac_03', 'mudfrac_04',         ...
        'mudfrac_05', 'mudfrac_06', 'mudfrac_07', 'mudfrac_08'}
    class = textscan(Vname, 'mudfrac_ %d');
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'Nbed';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['cohesive sediment fraction, ',        ...
                                 'size class ',                         ...
                                 num2str(class{:}, '%2.2i')];
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'ocean_time';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho Nbed ocean_time';
    else
      V.Attributes(3).Value   = 'x_rho y_rho Nbed ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case {'mudmass_01', 'mudmass_02', 'mudmass_03', 'mudmass_04',         ...
        'mudmass_05', 'mudmass_06', 'mudmass_07', 'mudmass_08'}
    class = textscan(Vname, 'mudmass_ %d');
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'Nbed';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['cohesive sediment mass, ',            ...
                                 'size class ',                         ...
                                 num2str(class{:}, '%2.2i')];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'kilogram meter-2';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho Nbed ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho Nbed ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case 'ripple_height'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'ocean_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'bottom ripple height';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case 'ripple_length'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'ocean_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'bottom ripple length';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
    
  case {'sand_01', 'sand_02', 'sand_03', 'sand_04',                     ...
        'sand_05', 'sand_06', 'sand_07', 'sand_08'}
    class = textscan(Vname, 'sand_ %d');
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['suspended noncohesive sediment, ',    ...
                                 'size class ',                         ...
                                 num2str(class{:}, '%2.2i')];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'kilogram meter-3';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho s_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho s_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case {'sandfrac_01', 'sandfrac_02', 'sandfrac_03', 'sandfrac_04',     ...
        'sandfrac_05', 'sandfrac_06', 'sandfrac_07', 'sandfrac_08'}
    class = textscan(Vname, 'sandfrac_ %d');
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'Nbed';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['noncohesive sediment fraction, ',     ...
                                 'size class ',                         ...
                                 num2str(class{:}, '%2.2i')];
    V.Attributes(2).Name      = 'time';
    V.Attributes(2).Value     = 'ocean_time';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho Nbed ocean_time';
    else
      V.Attributes(3).Value   = 'x_rho y_rho Nbed ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case {'sandmass_01', 'sandmass_02', 'sandmass_03', 'sandmass_04',     ...
        'sandmass_05', 'sandmass_06', 'sandmass_07', 'sandmass_08'}
    class = textscan(Vname, 'sandmass_ %d');
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'Nbed';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Dimensions(4).Name      = 'ocean_time';
    V.Dimensions(4).Length    = [];
    V.Dimensions(4).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = ['noncohesive sediment mass, ',         ...
                                 'size class ',                         ...
                                 num2str(class{:}, '%2.2i')];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'kilogram meter-2';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho Nbed ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho Nbed ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

  case 'settling_vel'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 'ocean_time';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = Unlimited;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'sediment median grain settling velocity';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'meter second-1';
    V.Attributes(3).Name      = 'time';
    V.Attributes(3).Value     = 'ocean_time';
    V.Attributes(4).Name      = 'coordinates';
    if (spherical),
      V.Attributes(4).Value   = 'lon_rho lat_rho ocean_time';
    else
      V.Attributes(4).Value   = 'x_rho y_rho ocean_time';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
    
%--------------------------------------------------------------------------
%  Inverse nudging time scales.
%--------------------------------------------------------------------------

  case 'M2_NudgeCoef'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = '2D momentum inverse nudging coefficients';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day-1';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho';
    else
      V.Attributes(3).Value   = 'x_rho y_rho';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
 
 case 'M3_NudgeCoef'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = '3D momentum inverse nudging coefficients';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day-1';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho s_rho';
    else
      V.Attributes(3).Value   = 'x_rho y_rho s_rho';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
 
 case 'tracer_NudgeCoef'
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = 'generic tracer inverse nudging coefficients';
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day-1';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho s_rho';
    else
      V.Attributes(3).Value   = 'x_rho y_rho s_rho';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);
 
 case {'temp_NudgeCoef',    'salt_NudgeCoef',                          ...
        'DON_NudgeCoef',     'iron_NudgeCoef',    'NO3_NudgeCoef',      ...
        'opal_NudgeCoef',    'PON_NudgeCoef',     'SiOH4_NudgeCoef',    ...
        'TIC_NudgeCoef',                                                ...    
        'mud_01_NudgeCoef',  'mud_02_NudgeCoef',  'mud_03_NudgeCoef',   ...
        'mud_04_NudgeCoef',  'mud_05_NudgeCoef',  'mud_06_NudgeCoef',   ...
        'mud_07_NudgeCoef',  'mud_08_NudgeCoef',                        ...
        'sand_01_NudgeCoef', 'sand_02_NudgeCoef', 'sand_03_NudgeCoef',  ...
        'sand_04_NudgeCoef', 'sand_05_NudgeCoef', 'sand_06_NudgeCoef',  ...
        'sand_07_NudgeCoef', 'sand_08_NudgeCoef'}
    field = Vname(1:strfind(Vname, '_NudgeCoef')-1);
    V.Name                    = Vname;
    V.Dimensions(1).Name      = 'xi_rho';
    V.Dimensions(1).Length    = [];
    V.Dimensions(1).Unlimited = false;
    V.Dimensions(2).Name      = 'eta_rho';
    V.Dimensions(2).Length    = [];
    V.Dimensions(2).Unlimited = false;
    V.Dimensions(3).Name      = 's_rho';
    V.Dimensions(3).Length    = [];
    V.Dimensions(3).Unlimited = false;
    V.Size                    = [];
    V.Attributes(1).Name      = 'long_name';
    V.Attributes(1).Value     = [field, ' inverse nudging coefficients'];
    V.Attributes(2).Name      = 'units';
    V.Attributes(2).Value     = 'day-1';
    V.Attributes(3).Name      = 'coordinates';
    if (spherical),
      V.Attributes(3).Value   = 'lon_rho lat_rho s_rho';
    else
      V.Attributes(3).Value   = 'x_rho y_rho s_rho';
    end
    V.Cgridtype.Name          = 'density point';
    V.Cgridtype.Value         = 1;
    V.Datatype                = Datatype;
    V.ncType                  = nc_constant(nctype);

%--------------------------------------------------------------------------
%  Requested variable not found.
%--------------------------------------------------------------------------

 otherwise
    disp(blanks(1));
    error(['ROMS_METADATA: unable to find metadata for variable: ',Vname]);
end

return
