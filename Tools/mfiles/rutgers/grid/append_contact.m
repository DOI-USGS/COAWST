function F=append_contact(Gnames, Cname, varargin)

%
% APPEND_CONTACT:  Appends nesting contact points data to ROMS grid
%
% F=append_contact(Gnames, Cname, Lplot);
%
% This function appends nesting contact point data to ROMS regular
% grid data. The resulting grid arrays are expanded to contain the
% data in the contact region.
%
% On Input:
%
%    Gnames        ROMS Grid NetCDF files/URL name (cell array)
%              or, an existing ROMS grid structure (struct array)
%
%    Cnames        ROMS nested grid connectivity NetCDF file (string)
%              or, an existing nesting connectivity structure (struct)
%
%    Lplot         Switch to plot various Contact Points figures
%                  (OPTIONAL, default false)
%
% On Output:
%
%    F             Structure containing expanded grid plus contact
%                  point data (struct array). The expanded grid coordinates
%                  depend on Cartesian or spherical grids.  The physical
%                  grid perimeter is also provided to aid in plotting.
%
%                  F(ng).x_rho(:,:)        or    F(ng).lon_rho(:,:) 
%                  F(ng).y_rho(:,:)        or    F(ng).lat_rho(:,:)
%                  F(ng).x_u(:,:)          or    F(ng).lon_u(:,:)
%                  F(ng).y_u(:,:)          or    F(ng).lat_u(:,:)
%                  F(ng).x_v(:,:)          or    F(ng).lon_v(:,:)
%                  F(ng).y_v(:,:)          or    F(ng).lat_v(:,:)
%                  F(ng).x_perimeter(:)    or    F(ng).lon_perimeter(:)
%                  F(ng).y_perimeter(:)    or    F(ng).lat_perimeter(:)
%
%                  F(ng).mask_rho(:,:)
%                  F(ng).mask_u(:,:)
%                  F(ng).mask_v(:,:)
%                  F(ng).h(:,:)
%                  F(ng).f(:,:)
%                  F(ng).pm(:,:)
%                  F(ng).pn(:,:)
%                  F(ng).angle(:,:)
%                  F(ng).dndx(:,:)
%                  F(ng).dmde(:,:)
%                  F(ng).area(:,:)
%                  F(ng).volume(:,:)
%
% Example:
%
%  Gnames = {'lake_jersey_grd_a.nc',  ...
%            'lake_jersey_grd_c.nc',  ...
%            'lake_jersey_grd_d.nc',  ...
%            'lake_jersey_grd_e.nc'};
%  Cname = 'lake_jersey_ngc_4g_acde.nc';
%
%  F = append_contact(Gnames, Cname);
%
%  ng = 2;
%  pcolor(F(ng).x_rho, F(ng).y_rho, 1./F(ng).pm); shading faceted;
%  hold on; plot(F(ng).x_perimeter, F(ng).y_perimeter, 'r-'); hold off
%

% svn $Id: append_contact.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

Lplot = false;

switch numel(varargin)
  case 1
    Lplot = varargin{1};
end

%--------------------------------------------------------------------------
% Set nested grids structure.
%--------------------------------------------------------------------------

if (~isstruct(Gnames)),
  G = grids_structure(Gnames);
else
  G = Gnames;
end

Ngrids = length(G);

%--------------------------------------------------------------------------
% Set nested grids connectivity structure.
%--------------------------------------------------------------------------

if (~isstruct(Cname)),
  S = read_contact(Cname);
else
  S = Cname;
end

Ncontact = S.Ncontact;

%--------------------------------------------------------------------------
% Build grid plust contact points data.
%--------------------------------------------------------------------------

% Allocate output structure.

if (G(1).spherical),
  F(1:Ngrids) = struct('lon_rho'      , [], 'lat_rho'      , [],        ...
                       'lon_u'        , [], 'lat_u'        , [],        ...
                       'lon_v'        , [], 'lat_v'        , [],        ...
                       'lon_perimeter', [], 'lat_perimeter', [],        ...
                       'mask_rho'     , [], 'mask_u'       , [],        ...
                       'mask_v'       , [],                             ...
                       'h'            , [], 'f'            , [],        ...
                       'pm'           , [], 'pn'           , [],        ...
                       'angle'        , [],                             ...
                       'dndx'         , [], 'dmde'         , [],        ...
                       'area'         , [], 'volume'       , []);
else
  F(1:Ngrids) = struct('x_rho'        , [], 'y_rho'        , [],        ...
                       'x_u'          , [], 'y_u'          , [],        ...
                       'x_v'          , [], 'y_v'          , [],        ...
                       'x_perimeter'  , [], 'y_perimeter'  , [],        ...
                       'mask_rho'     , [], 'mask_u'       , [],        ...
                       'mask_v'       , [],                             ...
                       'h'            , [], 'f'            , [],        ...
                       'pm'           , [], 'pn'           , [],        ...
                       'angle'        , [],                             ...
                       'dndx'         , [], 'dmde'         , [],        ...
                       'area'         , [], 'volume'       , []);
end  

% Fill grid information.

for ng=1:Ngrids,

  [ImR,JmR] = size(G(ng).mask_rho);
  [ImU,JmU] = size(G(ng).mask_u);
  [ImV,JmV] = size(G(ng).mask_v);

  IstrP = 1;   IendP = ImU;
  JstrP = 1;   JendP = JmV;

  if (ng == 1),
    if (G(ng).spherical),
      F(ng).lon_rho = G(ng).lon_rho;
      F(ng).lat_rho = G(ng).lat_rho;
      F(ng).lon_u   = G(ng).lon_u;
      F(ng).lat_u   = G(ng).lat_u;
      F(ng).lon_v   = G(ng).lon_v;
      F(ng).lat_v   = G(ng).lat_v;
    
      Xbox = [squeeze(G(ng).lon_psi(IstrP:IendP,JstrP));                ...
              squeeze(G(ng).lon_psi(IendP,JstrP+1:JendP))';             ...
              squeeze(flipud(G(ng).lon_psi(IstrP:IendP-1,JendP)));      ...
              squeeze(fliplr(G(ng).lon_psi(IstrP,JstrP:JendP-1)))'];

      Ybox = [squeeze(G(ng).lat_psi(IstrP:IendP,JstrP));                ...
              squeeze(G(ng).lat_psi(IendP,JstrP+1:JendP))';             ...
              squeeze(flipud(G(ng).lat_psi(IstrP:IendP-1,JendP)));      ...
              squeeze(fliplr(G(ng).lat_psi(IstrP,JstrP:JendP-1)))'];  
    
      F(ng).lon_perimeter = Xbox;
      F(ng).lat_perimeter = Ybox;
    else
      F(ng).x_rho = G(ng).x_rho;
      F(ng).y_rho = G(ng).y_rho;
      F(ng).x_u   = G(ng).x_u;
      F(ng).y_u   = G(ng).y_u;
      F(ng).x_v   = G(ng).x_v;
      F(ng).y_v   = G(ng).y_v;
    
      Xbox = [squeeze(G(ng).x_psi(IstrP:IendP,JstrP));                  ...
              squeeze(G(ng).x_psi(IendP,JstrP+1:JendP))';               ...
              squeeze(flipud(G(ng).x_psi(IstrP:IendP-1,JendP)));        ...
              squeeze(fliplr(G(ng).x_psi(IstrP,JstrP:JendP-1)))'];

      Ybox = [squeeze(G(ng).y_psi(IstrP:IendP,JstrP));                  ...
              squeeze(G(ng).y_psi(IendP,JstrP+1:JendP))';               ...
              squeeze(flipud(G(ng).y_psi(IstrP:IendP-1,JendP)));        ...
              squeeze(fliplr(G(ng).y_psi(IstrP,JstrP:JendP-1)))'];

      F(ng).x_perimeter = Xbox;
      F(ng).y_perimeter = Ybox;
    end
    F(ng).mask_rho = G(ng).mask_rho;
    F(ng).mask_u   = G(ng).mask_u;
    F(ng).mask_v   = G(ng).mask_v;
  
    F(ng).h     = G(ng).h;
    F(ng).f     = G(ng).f;
    F(ng).pm    = G(ng).pm;
    F(ng).pn    = G(ng).pn;
    F(ng).angle = G(ng).angle;
  
    if (~isempty(G(ng).dndx))
      F(ng).dndx  = G(ng).dndx;
    end
    
    if (~isempty(G(ng).dmde))
      F(ng).dmde  = G(ng).dmde;
    end  
  
  else
   
    if (G(ng).spherical),
      F(ng).lon_rho = nan(ImR+5,JmR+5);
      F(ng).lat_rho = nan(ImR+5,JmR+5);
      F(ng).lon_rho(4:ImR+3,4:JmR+3) = G(ng).lon_rho;
      F(ng).lat_rho(4:ImR+3,4:JmR+3) = G(ng).lat_rho;

      F(ng).lon_u = nan(ImU+5,JmU+5);
      F(ng).lat_u = nan(ImU+5,JmU+5);
      F(ng).lon_u(4:ImU+3,4:JmU+3) = G(ng).lon_u;
      F(ng).lat_u(4:ImU+3,4:JmU+3) = G(ng).lat_u;

      F(ng).lon_v = nan(ImV+5,JmV+5);
      F(ng).lat_v = nan(ImV+5,JmV+5);
      F(ng).lon_v(4:ImV+3,4:JmV+3) = G(ng).lon_v;
      F(ng).lat_v(4:ImV+3,4:JmV+3) = G(ng).lat_v;
    
      Xbox = [squeeze(G(ng).lon_psi(IstrP:IendP,JstrP));                ...
              squeeze(G(ng).lon_psi(IendP,JstrP+1:JendP))';             ...
              squeeze(flipud(G(ng).lon_psi(IstrP:IendP-1,JendP)));      ...
              squeeze(fliplr(G(ng).lon_psi(IstrP,JstrP:JendP-1)))'];

      Ybox = [squeeze(G(ng).lat_psi(IstrP:IendP,JstrP));                ...
              squeeze(G(ng).lat_psi(IendP,JstrP+1:JendP))';             ...
              squeeze(flipud(G(ng).lat_psi(IstrP:IendP-1,JendP)));      ...
              squeeze(fliplr(G(ng).lat_psi(IstrP,JstrP:JendP-1)))'];  
   
      F(ng).lon_perimeter = Xbox;
      F(ng).lat_perimeter = Ybox;
    else
      F(ng).x_rho = nan(ImR+5,JmR+5);
      F(ng).y_rho = nan(ImR+5,JmR+5);
      F(ng).x_rho(4:ImR+3,4:JmR+3) = G(ng).x_rho;
      F(ng).y_rho(4:ImR+3,4:JmR+3) = G(ng).y_rho;

      F(ng).x_u = nan(ImU+5,JmU+5);
      F(ng).y_u = nan(ImU+5,JmU+5);
      F(ng).x_u(4:ImU+3,4:JmU+3) = G(ng).x_u;
      F(ng).y_u(4:ImU+3,4:JmU+3) = G(ng).y_u;

      F(ng).x_v = nan(ImV+5,JmV+5);
      F(ng).y_v = nan(ImV+5,JmV+5);
      F(ng).x_v(4:ImV+3,4:JmV+3) = G(ng).x_v;
      F(ng).y_v(4:ImV+3,4:JmV+3) = G(ng).y_v;
    
      Xbox = [squeeze(G(ng).x_psi(IstrP:IendP,JstrP));                  ...
              squeeze(G(ng).x_psi(IendP,JstrP+1:JendP))';               ...
              squeeze(flipud(G(ng).x_psi(IstrP:IendP-1,JendP)));        ...
              squeeze(fliplr(G(ng).x_psi(IstrP,JstrP:JendP-1)))'];

      Ybox = [squeeze(G(ng).y_psi(IstrP:IendP,JstrP));                  ...
              squeeze(G(ng).y_psi(IendP,JstrP+1:JendP))';               ...
              squeeze(flipud(G(ng).y_psi(IstrP:IendP-1,JendP)));        ...
              squeeze(fliplr(G(ng).y_psi(IstrP,JstrP:JendP-1)))'];

      F(ng).x_perimeter = Xbox;
      F(ng).y_perimeter = Ybox;
    end
    F(ng).mask_rho = nan(ImR+5,JmR+5);
    F(ng).mask_u   = nan(ImU+5,JmU+5);
    F(ng).mask_v   = nan(ImV+5,JmV+5);

    F(ng).mask_rho(4:ImR+3,4:JmR+3) = G(ng).mask_rho;
    F(ng).mask_u(4:ImU+3,4:JmU+3)   = G(ng).mask_u;
    F(ng).mask_v(4:ImV+3,4:JmV+3)   = G(ng).mask_v;

    F(ng).h     = nan(ImR+5,JmR+5);
    F(ng).f     = nan(ImR+5,JmR+5);
    F(ng).pm    = nan(ImR+5,JmR+5);
    F(ng).pn    = nan(ImR+5,JmR+5);
    F(ng).angle = nan(ImR+5,JmR+5);

    F(ng).h(4:ImR+3,4:JmR+3)     = G(ng).h;
    F(ng).f(4:ImR+3,4:JmR+3)     = G(ng).f;
    F(ng).pm(4:ImR+3,4:JmR+3)    = G(ng).pm;
    F(ng).pn(4:ImR+3,4:JmR+3)    = G(ng).pn;
    F(ng).angle(4:ImR+3,4:JmR+3) = G(ng).angle;

    if (~isempty(G(ng).dndx))
      F(ng).dndx  = nan(ImR+5,JmR+5);
      F(ng).dndx(4:ImR+3,4:JmR+3)  = G(ng).dndx;
    else
      F(ng).dndx  = zeros(ImR+5,JmR+5);
    end
    
    if (~isempty(G(ng).dmde))
      F(ng).dmde  = nan(ImR+5,JmR+5);
      F(ng).dmde(4:ImR+3,4:JmR+3)  = G(ng).dmde;
    else
      F(ng).dmde  = zeros(ImR+5,JmR+5);
    end
  
  end

end

% Fill contact points values.

do_grid = true([1 Ngrids]);

for cr=1:Ncontact,
  rg = S.contact(cr).receiver_grid;

  if (do_grid(rg) && rg > 1),

    do_grid(rg) = false;

    NpointsR = length(S.contact(cr).point.Irg_rho);
    for n=1:NpointsR,
      Ir = S.contact(cr).point.Irg_rho(n)+4;
      Jr = S.contact(cr).point.Jrg_rho(n)+4;
      if (G(rg).spherical),
        F(rg).lon_rho(Ir,Jr) = S.contact(cr).point.Xrg_rho(n);
        F(rg).lat_rho(Ir,Jr) = S.contact(cr).point.Yrg_rho(n);
      else
        F(rg).x_rho(Ir,Jr) = S.contact(cr).point.Xrg_rho(n);
        F(rg).y_rho(Ir,Jr) = S.contact(cr).point.Yrg_rho(n);
      end
      F(rg).mask_rho(Ir,Jr) = S.contact(cr).point.mask_rho(n);

      F(rg).h(Ir,Jr)     = S.contact(cr).point.h(n);
      F(rg).f(Ir,Jr)     = S.contact(cr).point.f(n);
      F(rg).pm(Ir,Jr)    = S.contact(cr).point.pm(n);
      F(rg).pn(Ir,Jr)    = S.contact(cr).point.pn(n);
      F(rg).angle(Ir,Jr) = S.contact(cr).point.angle(n);
    
      if (~isempty(G(rg).dndx))
        F(rg).dndx(Ir,Jr)  = S.contact(cr).point.dndx(n);
      end
      
      if (~isempty(G(rg).dmde))
        F(rg).dmde(Ir,Jr)  = S.contact(cr).point.dmde(n);
      end      
    
    end

    NpointsU = length(S.contact(cr).point.Irg_u);
    for n=1:NpointsU,
      Iu = S.contact(cr).point.Irg_u(n)+3;
      Ju = S.contact(cr).point.Jrg_u(n)+4;
      if (G(rg).spherical),
        F(rg).lon_u(Iu,Ju) = S.contact(cr).point.Xrg_u(n);
        F(rg).lat_u(Iu,Ju) = S.contact(cr).point.Yrg_u(n);
      else
        F(rg).x_u(Iu,Ju) = S.contact(cr).point.Xrg_u(n);
        F(rg).y_u(Iu,Ju) = S.contact(cr).point.Yrg_u(n);
      end
      F(rg).mask_u(Iu,Ju) = S.contact(cr).point.mask_u(n);
    end

    NpointsV = length(S.contact(cr).point.Irg_v);
    for n=1:NpointsV,
      Iv = S.contact(cr).point.Irg_v(n)+4;
      Jv = S.contact(cr).point.Jrg_v(n)+3;
      F(rg).mask_v(Iv,Jv) = S.contact(cr).point.mask_v(n);
      if (G(rg).spherical),
        F(rg).lon_v(Iv,Jv) = S.contact(cr).point.Xrg_v(n);
        F(rg).lat_v(Iv,Jv) = S.contact(cr).point.Yrg_v(n);
      else
        F(rg).x_v(Iv,Jv) = S.contact(cr).point.Xrg_v(n);
        F(rg).y_v(Iv,Jv) = S.contact(cr).point.Yrg_v(n);
      end
    end

  end

end

% Compute area and volume.

for ng=1:Ngrids,
  dx = 1./F(ng).pm;
  dy = 1./F(ng).pn;
  F(ng).area = dx .* dy;
  F(ng).volume = F(ng).area .* F(ng).h;
end

%--------------------------------------------------------------------------
%  Plot various grid fields.
%--------------------------------------------------------------------------

if (Lplot),
  cst = 'm-';
  box = 'r-';

  Mversion = version('-release');
  Vyear    = sscanf(Mversion, '%i');

  if (strcmp(Mversion,'2014b') || Vyear > 2014),
    h = gcf;
    nf = h.Number;
  else
    ag = findobj;                                   % all graphical objects
    nf = max(ag(ag == fix(ag)));                    % number of figures
  end
    
  for ng=1:Ngrids,
    if (G(ng).spherical)
      Xr = F(ng).lon_rho;
      Yr = F(ng).lat_rho;
      Xperimeter = F(ng).lon_perimeter;
      Yperimeter = F(ng).lat_perimeter;
      if (isfield(G(ng),'lon_coast') && isfield(G(ng),'lat_coast')),
        Clon = G(ng).lon_coast;
        Clat = G(ng).lat_coast;
        got.coast = true;
      end
    else
      Xr = 0.001 .* F(ng).x_rho;                    % km
      Yr = 0.001 .* F(ng).y_rho;
      Xperimeter = 0.001 .* F(ng).x_perimeter;
      Yperimeter = 0.001 .* F(ng).y_perimeter;
      got.coast = true;
    end
  
    ic = nf;

    if (ng == 1),
      Xmin=min(Xr(:));
      Xmax=max(Xr(:));
      Ymin=min(Yr(:));
      Ymax=max(Yr(:));
      
      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Bathymetry: Interior plus Contact Points');

      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Coriolis: Interior plus Contact Points');
      
      if (length(unique(F(1).angle(:))) ~= 1),
        ic = ic + 1;
        figure(ic);
        plot(Xperimeter, Yperimeter, box);
        axis([Xmin Xmax Ymin Ymax]); grid on;
        hold on;
        title('Curvilinear Angle: Interior plus Contact Points');
      end

      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Inverse Grid X-Spacing: Interior plus Contact Points');

      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Inverse Grid Y-Spacing: Interior plus Contact Points');

      if (length(unique(F(1).dndx(:))) ~= 1),
        ic = ic + 1;
        figure(ic);
        plot(Xperimeter, Yperimeter, box);
        axis([Xmin Xmax Ymin Ymax]); grid on;
        hold on;
        title('Curvilinear Metric DNDX: Interior plus Contact Points');
      end

      if (length(unique(F(1).dmde(:))) ~= 1),
        ic = ic + 1;
        figure(ic);
        plot(Xperimeter, Yperimeter, box);
        axis([Xmin Xmax Ymin Ymax]); grid on;
        hold on;
        title('Curvilinear Metric DMDE: Interior plus Contact Points');
      end
      
      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Grid Area: Interior plus Contact Points');
      
      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Grid Volume: Interior plus Contact Points');
    else
      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xr, Yr, F(ng).h); shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      if (got.coast && ng == Ngrids)
        plot(Clon, Clat, cst);
      end

      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xr, Yr, F(ng).f); shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      if (got.coast && ng == Ngrids)
        plot(Clon, Clat, cst);
      end

      if (length(unique(F(1).angle(:))) ~= 1),
        ic = ic + 1;
        figure(ic);
        pcolorjw (Xr, Yr, F(ng).angle); shading flat; colorbar;
        plot(Xperimeter, Yperimeter, box);
        if (got.coast && ng == Ngrids)
          plot(Clon, Clat, cst);
        end
      end

      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xr, Yr, F(ng).pm); shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      if (got.coast && ng == Ngrids)
        plot(Clon, Clat, cst);
      end

      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xr, Yr, F(ng).pn); shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      if (got.coast && ng == Ngrids)
        plot(Clon, Clat, cst);
      end
      
      if (length(unique(F(1).dndx(:))) ~= 1),
        ic = ic + 1;
        figure(ic);
        pcolorjw (Xr, Yr, F(ng).dndx); shading flat; colorbar;
        plot(Xperimeter, Yperimeter, box);
        if (got.coast && ng == Ngrids)
          plot(Clon, Clat, cst);
        end
      end

      if (length(unique(F(1).dmde(:))) ~= 1),
        ic = ic + 1;
        figure(ic);
        pcolorjw (Xr, Yr, F(ng).dmde); shading flat; colorbar;
        plot(Xperimeter, Yperimeter, box);
        if (got.coast && ng == Ngrids)
          plot(Clon, Clat, cst);
        end
      end
      
      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xr, Yr, F(ng).area); shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      if (got.coast && ng == Ngrids)
        plot(Clon, Clat, cst);
      end

      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xr, Yr, F(ng).volume); shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      if (got.coast && ng == Ngrids)
        plot(Clon, Clat, cst);
      end
    end
  end
end

return
