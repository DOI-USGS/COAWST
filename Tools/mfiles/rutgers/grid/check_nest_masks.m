function [S, G] = check_nest_masks(Gnames, Cname, Lplot)

%
% CHECK_NEST_MASKS:  Checks land/sea masks in nesting applications
%
% [S, G] = check_nest_masks(Gnames, Cname, Lplot);
%
% This function can be used to check land/sea masking over nesting
% contact points. In particular, it can be used to examine the nesting
% grids connectivity when there are land/sea mask features in the contact
% regions.
%
% The plotting section below is provided as a guideline.  Users may
% need to adapt it for their particular application.
%
% On Input:
%
%    Gnames        ROMS Grid/History NetCDF file/URL names containing
%                    all grid variables (cell array)
%
%    Cname         Contact Points NetCDF file name (string), generated
%                    from "contact.m".
%
%    Lplot         Switch to plot various land/sea maskinf figures
%                    (default true)
%
% On Output:
%
%    S             Nested grids Contact Points structure (struct array)
%
%    G             Nested grids structure (1 x Ngrids struct array).
%                    Additional fields are added to the structure to
%                    include regular grid plus contact points:
%
%                    G(ng).nest_lon_rho
%                    G(ng).nest_lat_rho
%                    G(ng).nest_lon_u
%                    G(ng).nest_lat_u
%                    G(ng).nest_lon_v
%
%                    G(ng).nest_mask_rho
%                    G(ng).nest_mask_u
%                    G(ng).nest_mask_v
%
% Example:
%
% [S,G] = check_nest_masks({'my_grid_coarser.nc', 'my_grid_finer.nc',}, ...
%                          'my_contact.nc', true);
%
% Calls to External Functions:
%
%    grids_structure     Gets nested grid Structure, G(ng) 
%    read_contact        Gets nesting grid connectivity structure, S
%

% svn $Id: check_nest_masks.m 732 2014-04-18 18:24:32Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

if (nargin < 2),
  Lplot = true;
end

%--------------------------------------------------------------------------
% Get nested grids structure.
%--------------------------------------------------------------------------

G = grids_structure(Gnames);

Ngrids = length(G);

%--------------------------------------------------------------------------
% Read and set nesting grids contact points structure.
%--------------------------------------------------------------------------

S = read_contact(Cname);

Ncontact = S.Ncontact;

%--------------------------------------------------------------------------
% Build mask with contact points.
%--------------------------------------------------------------------------

for ng=1:Ngrids,
  if (G(ng).spherical),
    G(ng).nest_lon_rho = [];
    G(ng).nest_lat_rho = [];
    G(ng).nest_lon_u   = [];
    G(ng).nest_lat_u   = [];
    G(ng).nest_lon_v   = [];
    G(ng).nest_lat_v   = [];
  else
    G(ng).nest_x_rho = [];
    G(ng).nest_y_rho = [];
    G(ng).nest_x_u   = [];
    G(ng).nest_y_u   = [];
    G(ng).nest_x_v   = [];
    G(ng).nest_y_v   = [];
  end
  G(ng).nest_mask_rho = [];
  G(ng).nest_mask_u   = [];
  G(ng).nest_mask_v   = [];
end

for cr=1:Ncontact,
  dg = S.contact(cr).donor_grid;
  rg = S.contact(cr).receiver_grid;

  if (rg > 1),

% Mask at RHO-points.

    [Im,Jm] = size(G(rg).mask_rho);

    G(rg).nest_mask_rho = nan(Im+5,Jm+5);
    G(rg).nest_mask_rho(4:Im+3,4:Jm+3) = G(rg).mask_rho;

    if (G(rg).spherical),
      G(rg).nest_lon_rho = nan(Im+5,Jm+5);
      G(rg).nest_lat_rho = nan(Im+5,Jm+5);
      G(rg).nest_lon_rho(4:Im+3,4:Jm+3) = G(rg).lon_rho;
      G(rg).nest_lat_rho(4:Im+3,4:Jm+3) = G(rg).lat_rho;
    else
      G(rg).nest_x_rho = nan(Im+5,Jm+5);
      G(rg).nest_y_rho = nan(Im+5,Jm+5);
      G(rg).nest_x_rho(4:Im+3,4:Jm+3) = G(rg).x_rho;
      G(rg).nest_y_rho(4:Im+3,4:Jm+3) = G(rg).y_rho;   
    end
    
    Npoints = length(S.contact(cr).point.Irg_rho);
    for n=1:Npoints,
      Ir = S.contact(cr).point.Irg_rho(n)+4;
      Jr = S.contact(cr).point.Jrg_rho(n)+4;
      G(rg).nest_mask_rho(Ir,Jr) = S.contact(cr).point.mask_rho(n);
      if (G(ng).spherical),
        G(rg).nest_lon_rho(Ir,Jr) = S.contact(cr).point.Xrg_rho(n);
        G(rg).nest_lat_rho(Ir,Jr) = S.contact(cr).point.Yrg_rho(n);     
      else
        G(rg).x_rho(Ir,Jr) = S.contact(cr).point.Xrg_rho(n);
        G(rg).y_rho(Ir,Jr) = S.contact(cr).point.Yrg_rho(n);      
      end
    end

% Mask at U-points.

    [Im,Jm] = size(G(rg).mask_u);

    G(rg).nest_mask_u = nan(Im+5,Jm+5);
    G(rg).nest_mask_u(4:Im+3,4:Jm+3) = G(rg).mask_u;

    if (G(rg).spherical),
      G(rg).nest_lon_u = nan(Im+5,Jm+5);
      G(rg).nest_lat_u = nan(Im+5,Jm+5);
      G(rg).nest_lon_u(4:Im+3,4:Jm+3) = G(rg).lon_u;
      G(rg).nest_lat_u(4:Im+3,4:Jm+3) = G(rg).lat_u;
    else
      G(rg).nest_x_u = nan(Im+5,Jm+5);
      G(rg).nest_y_u = nan(Im+5,Jm+5);
      G(rg).nest_x_u(4:Im+3,4:Jm+3) = G(rg).x_u;
      G(rg).nest_y_u(4:Im+3,4:Jm+3) = G(rg).y_u;   
    end
    
    Npoints = length(S.contact(cr).point.Irg_u);
    for n=1:Npoints,
      Iu = S.contact(cr).point.Irg_u(n)+3;
      Ju = S.contact(cr).point.Jrg_u(n)+4;
      G(rg).nest_mask_u(Iu,Ju) = S.contact(cr).point.mask_u(n);
      if (G(ng).spherical),
        G(rg).nest_lon_u(Iu,Ju) = S.contact(cr).point.Xrg_u(n);
        G(rg).nest_lat_u(Iu,Ju) = S.contact(cr).point.Yrg_u(n);     
      else
        G(rg).x_u(Iu,Ju) = S.contact(cr).point.Xrg_u(n);
        G(rg).y_u(Iu,Ju) = S.contact(cr).point.Yrg_u(n);      
      end
    end

% Mask at V-points.

    [Im,Jm] = size(G(rg).mask_v);

    G(rg).nest_mask_v = nan(Im+5,Jm+5);
    G(rg).nest_mask_v(4:Im+3,4:Jm+3) = G(rg).mask_v;

    if (G(rg).spherical),
      G(rg).nest_lon_v = nan(Im+5,Jm+5);
      G(rg).nest_lat_v = nan(Im+5,Jm+5);
      G(rg).nest_lon_v(4:Im+3,4:Jm+3) = G(rg).lon_v;
      G(rg).nest_lat_v(4:Im+3,4:Jm+3) = G(rg).lat_v;
    else
      G(rg).nest_x_v = nan(Im+5,Jm+5);
      G(rg).nest_y_v = nan(Im+5,Jm+5);
      G(rg).nest_x_v(4:Im+3,4:Jm+3) = G(rg).x_v;
      G(rg).nest_y_v(4:Im+3,4:Jm+3) = G(rg).y_v;   
    end
    
    Npoints = length(S.contact(cr).point.Irg_v);
    for n=1:Npoints,
      Iv = S.contact(cr).point.Irg_v(n)+4;
      Jv = S.contact(cr).point.Jrg_v(n)+3;
      G(rg).nest_mask_v(Iv,Jv) = S.contact(cr).point.mask_v(n);
      if (G(ng).spherical),
        G(rg).nest_lon_v(Iv,Jv) = S.contact(cr).point.Xrg_v(n);
        G(rg).nest_lat_v(Iv,Jv) = S.contact(cr).point.Yrg_v(n);     
      else
        G(rg).x_v(Iv,Jv) = S.contact(cr).point.Xrg_v(n);
        G(rg).y_v(Iv,Jv) = S.contact(cr).point.Yrg_v(n);      
      end
    end

  end  

end

%--------------------------------------------------------------------------
% Plot land/sea masking fields.
%--------------------------------------------------------------------------

if (Lplot),
  
% RHO-points.

  for ng=1:Ngrids,
    if (ng == 1),
      figure;
      hr(ng) = pcolorjw(G(ng).lon_rho, G(ng).lat_rho, G(ng).mask_rho);
      title('Land/Sea Masking at RHO-points');
      shading faceted; colormap(gray);
      hold on;

      set(gcf, 'renderer', 'OpenGL');
      set(hr(ng), 'EdgeColor', [0.4 0.4 0.4]);
      alpha(hr(ng), 0.2);
    else
%     hr(ng) = pcolorjw(G(ng).lon_rho, G(ng).lat_rho, G(ng).mask_rho);
      hr(ng) = pcolorjw(G(ng).nest_lon_rho, G(ng).nest_lat_rho,         ...
                        G(ng).nest_mask_rho);
      shading faceted; colormap(gray);
      xlabel('Contact points are included')
      set(hr(ng), 'EdgeColor', [0.7 0.7 0.7]);
      alpha(hr(ng), 0.2);
    end
  end

  saveas(gcf, 'rho_masks', 'fig');

% U-points.

  for ng=1:Ngrids,
    if (ng == 1),
      figure;
      hu(ng) = pcolorjw(G(ng).lon_u, G(ng).lat_u, G(ng).mask_u);
      title('Land/Sea Masking at U-points');
      shading faceted; colormap(gray);
      hold on;

      set(gcf, 'renderer', 'OpenGL');
      set(hu(ng), 'EdgeColor', [0.4 0.4 0.4]);
      alpha(hu(ng), 0.2);
    else
%     hu(ng) = pcolorjw(G(ng).lon_u, G(ng).lat_u, G(ng).mask_u);
      hu(ng) = pcolorjw(G(ng).nest_lon_u, G(ng).nest_lat_u,             ...
                        G(ng).nest_mask_u);
      shading faceted; colormap(gray);
      xlabel('Contact points are included')
      set(hu(ng), 'EdgeColor', [0.7 0.7 0.7]);
      alpha(hu(ng), 0.2);
    end
  end

  saveas(gcf, 'u_masks', 'fig');

% V-points.

  for ng=1:Ngrids,
    if (ng == 1),
      figure;
      hv(ng) = pcolorjw(G(ng).lon_v, G(ng).lat_v, G(ng).mask_v);
      title('Land/Sea Masking at V-points');
      shading faceted; colormap(gray);
      hold on;

      set(gcf, 'renderer', 'OpenGL');
      set(hv(ng), 'EdgeColor', [0.4 0.4 0.4]);
      alpha(hv(ng), 0.2);
    else
%     hv(ng) = pcolorjw(G(ng).lon_v, G(ng).lat_v, G(ng).mask_v);
      hv(ng) = pcolorjw(G(ng).nest_lon_v, G(ng).nest_lat_v,             ...
                        G(ng).nest_mask_v);
      shading faceted; colormap(gray);
      xlabel('Contact points are included')
      set(hv(ng), 'EdgeColor', [0.7 0.7 0.7]);
      alpha(hv(ng), 0.2);
    end
  end

  saveas(gcf, 'v_masks', 'fig');

end

return