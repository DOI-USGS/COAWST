function G = flip_grid (Ginp, Gout)
  
%
% GRID_FLIP:  Creates a ROMS Grid NetCDF file with flipped coordinates
%
% G = flip_grid (Ginp, Gout)
%
% Given a ROMS Grid NetCDF file (Ginp), this function creates a new
% Grid NetCDF file (Gout) with the flipped dimensions, coordinates,
% and variables.
%
% On Input:
%
%    Ginp       Input ROMS Grid NetCDF file name (string)
%    Gout       Output flipped ROMS Grid NetCDF file name (string)
%
% On Output:
%
%    G          Flipped Grid structure (struct array)
%
% For example, the flipping is as follows:
%
%                 Output file          Input file   
%
%                 x_rho(:,:)     =     y_rho(:,:)'
%                 y_rho(:,:)     =     x_rho(:,:)'
%
%                 x_psi(:,:)     =     y_psi(:,:)'
%                 y_psi(:,:)     =     x_psi(:,:)'
%
%                 x_u(:,:)       =     y_v(:,:)'
%                 y_u(:,:)       =     x_v(:,:)'
%
%                 x_v(:,:)       =     y_u(:,:)'
%                 y_v(:,:)       =     x_u(:,:)'
%
%                 mask_rho(:,:)  =     mask_rho(:,:)'
%
%                 pm(:,:)        =     pn(:,:)'
%                 pn(:,:)        =     pm(:,:)'
%

% svn $Id: flip_grid.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
  
% Inquire about the contents of the NetCDF file.

I = nc_inq(Ginp);

% Get input grid structure.

G = get_roms_grid (Ginp);

% Flip dimension and overwrite the values in the inquire structure, I.

Lr = I.Dimensions(strcmp({I.Dimensions.Name},'xi_rho' )).Length;
Mr = I.Dimensions(strcmp({I.Dimensions.Name},'eta_rho')).Length;

L = Mr-1;
M = Lr-1;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_rho' )).Length = L+1;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_rho')).Length = M+1;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_psi' )).Length = L;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_psi')).Length = M;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_u'   )).Length = L;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_u'  )).Length = M+1;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_v'   )).Length = L+1;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_v'  )).Length = M;

%--------------------------------------------------------------------------
% Create subdomain NetCDF, use native Matlab interface.
%--------------------------------------------------------------------------

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

nc_create(Gout,mode,I);

%--------------------------------------------------------------------------
% Write out subdomain NetCDF data, use native Matlab interface.
%--------------------------------------------------------------------------

nvars = length(I.Variables);

for n=1:nvars,
  vname = char(I.Variables(n).Name);
  vsize = length(I.Variables(n).Size);
  
  if (isempty(I.Variables(n).Dimensions)),

    Vinp = ncread(Ginp, vname);                 % process scalar data
    if (~isempty(Vinp)),
      ncwrite(Gout, vname, Vinp);
    end
  
  else

    unlimited = [I.Variables(n).Dimensions.Unlimited];
    if (any(unlimited)),
      mrec = I.Variables(n).Dimensions(unlimited).Length;
    else
      mrec = 0;
    end

    if (vsize >= 2),
      switch (vname)
        case 'lon_rho'
          ncwrite(Gout, vname, G.lat_rho');
        case 'lat_rho'
          ncwrite(Gout, vname, G.lon_rho');
        case 'x_rho'
          ncwrite(Gout, vname, G.y_rho');
        case 'y_rho'
          ncwrite(Gout, vname, G.x_rho');
        case 'lon_psi'
          ncwrite(Gout, vname, G.lat_psi');
        case 'lat_psi'
          ncwrite(Gout, vname, G.lon_psi');
        case 'x_psi'
          ncwrite(Gout, vname, G.y_psi');
        case 'y_psi'
          ncwrite(Gout, vname, G.x_psi');
        case 'lon_u'
          ncwrite(Gout, vname, G.lat_v');
        case 'lat_u'
          ncwrite(Gout, vname, G.lon_v');
        case 'x_u'
          ncwrite(Gout, vname, G.y_v');
        case 'y_u'
          ncwrite(Gout, vname, G.x_v');
        case 'lon_v'
          ncwrite(Gout, vname, G.lat_u');
        case 'lat_v'
          ncwrite(Gout, vname, G.lon_u');
        case 'x_v'
          ncwrite(Gout, vname, G.y_u');
        case 'y_v'
          ncwrite(Gout, vname, G.x_u');
        case 'mask_rho'
          ncwrite(Gout, vname, G.mask_rho');
        case 'mask_psi'
          ncwrite(Gout, vname, G.mask_psi');
        case 'mask_u'
          ncwrite(Gout, vname, G.mask_v');
        case 'mask_v'
          ncwrite(Gout, vname, G.mask_u');
        case 'pm'
          ncwrite(Gout, vname, G.pn');
        case 'pn'
          ncwrite(Gout, vname, G.pm');
        case 'hraw'
          for m=1:mrec,
            Vinp = ncread(Ginp, vname, [1 1 m], [Inf Inf 1]);
            ncwrite(Gout, vname, Vinp');
          end
        otherwise
          Vinp = ncread(Ginp, vname);
          if (vsize == 2),
            Vinp = permute(Vinp, [2,1]);
          elseif (vsize == 3),
            Vinp = permute(Vinp, [2,1,3]);
          end
          ncwrite(Gout, vname, Vinp);
      end      
    else
      Vinp  = ncread(Ginp, vname);
      ncwrite(Gout, vname, Vinp);
    end
  end
end

% Recompute Land/sea masking.

if (isfield(G,'mask_rho')),
  [umask, vmask, pmask] = uvp_masks(G.mask_rho');

  ncwrite(Gout, 'mask_psi', pmask);
  ncwrite(Gout, 'mask_u'  , umask);
  ncwrite(Gout, 'mask_v'  , vmask);
end

% Flip domain dimensions.

xl = ncread(Ginp, 'xl');
el = ncread(Ginp, 'el');

ncwrite(Gout, 'xl', el);
ncwrite(Gout, 'el', xl);

%--------------------------------------------------------------------------
% Get new grid structure.
%--------------------------------------------------------------------------

G = get_roms_grid (Gout);

return
