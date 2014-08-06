function I = nc_check(Info)

%
% NC_CHECK:  Checks NetCDF information structure for compliance
%
% I = nc_check(Iold)
%
% This function checks the information structure returned from calls to
% "nc_inq" or native "ncinfo" for compliance and changes new variable
% types and attributes. For example, it updates the "spherical" switch
% to integer and fixes the land/sea masking attributes for compliance.
%
% Recall that this information structure is very convenient to create
% new NetCDF files having the dimensions and variable schema of and
% old NetCDF file.
%
% On Input:
%
%    Info        NetCDF file information structure (struct array)
%
% On Ouput:
%
%    I           Updated NetCDF file information structure (struct array)
%

% svn $Id: nc_check.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Note: A easy way to delete all attributes for a variable is to set its
% ====  value to empty:
%
%      I.Variables(index).Attributes = []
%
% or overwrite with an Attribute substructure, A:
%
%      I.Variables(index).Attributes = A
%

% Initialize output structure with input structure.
  
I = Info;

%--------------------------------------------------------------------------
% Check spherical variable and change to integer.  Nowadays, we avoid
% using a character.
%--------------------------------------------------------------------------

index = strcmp({I.Variables.Name}, 'spherical');

if (any(index) && strcmpi(I.Variables(index).Datatype, 'char')),

  I.Variables(index).Datatype   = 'int32';
  I.Variables(index).ncType     = nc_constant('nc_int');

  A(1).Name  = 'long_name';
  A(1).Value = 'grid type logical switch';

  A(2).Name  = 'flag_values';
  A(2).Value = [int32(0) int32(1)];
  
  A(3).Name  = 'flag_meanings';
  A(3).Value = 'Cartesian spherical';
  
  I.Variables(index).Attributes = A;
end

%--------------------------------------------------------------------------
% Check Land/Sea  masking attributes.
%--------------------------------------------------------------------------

Pindex = strcmp({I.Variables.Name}, 'mask_psi');
Rindex = strcmp({I.Variables.Name}, 'mask_rho');
Uindex = strcmp({I.Variables.Name}, 'mask_u');
Vindex = strcmp({I.Variables.Name}, 'mask_v');

A(1).Name  = 'long_name';
A(1).Value = 'mask on';

A(2).Name  = 'flag_values';
A(2).Value = [0 1];
  
A(3).Name  = 'flag_meanings';
A(3).Value = 'land water';
  
if (any(Pindex)),
  A(1).Value = 'mask on PSI-points';
  I.Variables(Pindex).Attributes = A;
end

if (any(Rindex)),
  A(1).Value = 'mask on RHO-points';
  I.Variables(Rindex).Attributes = A;
end

if (any(Uindex)),
  A(1).Value = 'mask on U-points';
  I.Variables(Uindex).Attributes = A;
end

if (any(Vindex)),
  A(1).Value = 'mask on V-points';
  I.Variables(Vindex).Attributes = A;
end

return
