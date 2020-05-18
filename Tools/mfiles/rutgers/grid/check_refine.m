function S = check_refine(Ginp, Cinp, varargin)

%
% CHECK_REFINE: Check a ROMS nested refinement grid structure
%
% S = check_refine(Ginp, Cinp, Lplot)
%
% Generates and plots several information fields to check a refinement
% nested configuration.  Several field are averaged in the fine-to-coarse
% sense to check conservation properties:
%
%       AVG[Fvalues(:,:)_i , i=1,rfactor^2] = Cvalues(Idg,Jdg)
%
% where (Idg,Jdg) are the indices of the coarse grid donor cell contatining
% the finer grid values.
%
% On Input:
%
%    Ginp          ROMS nested Grid NetCDF files/URL name (cell array)
%              or, an existing ROMS grid structure (struct array)
%
%    Cinp          ROMS nested grid connectivity NetCDF file (string)
%              or, an existing nesting connectivity structure (struct)
%
%    Lplot         Switch to plot various differences fields between
%                    coarser donor and averaged refined grid values
%                    (default false)
%
% On Output:
%
%    S             Interpolation information structure (struct array,
%                    OPTIONAL):
%
%                    S(ng).grid_name       grid NetCDF file name
%                    S(ng).spherical       spherical switch
%                    S(ng).refine_factor   grid refinement factor
%                    S(ng).parent_Imin     donor grid lower-left  I-index
%                    S(ng).parent_Imax     donor grid upper-right I-index
%                    S(ng).parent_Jmin     donor grid lower-left  J-index
%                    S(ng).parent_Jmax     donor grid upper-right J-index
%                    S(ng).x               Cartesian/spherical X-coordinate
%                    S(ng).y               Cartesian/spherical Y-coordinate
%                    S(ng).Idg_avg         averaged interior donor I-index
%                    S(ng).Jdg_avg         averaged interior donor J-index
%                    S(ng).dg_index        averaged interior donor 2D index
%                    S(ng).x_avg           averaged interior X-coordinate
%                    S(ng).y_avg           averaged interior Y-coordinate
%                    S(ng).x_perimeter     X-perimeter at PSI-points
%                    S(ng).y_perimeter     Y-perimeter at PSI-points
%                    S(ng).f               Coriolis parameter
%                    S(ng).f_avg           interior averaged "f"
%                    S(ng).f_diff          coarse minus averaged "f"
%                    S(ng).h               bathymetry
%                    S(ng).h_avg           interior averaged "h"
%                    S(ng).h_diff          coarse minus averaged "h"
%                    S(ng).angle           curvilinear angle
%                    S(ng).angle_avg       interior averaged "angle"
%                    S(ng).angle_diff      coarse minus averaged "angle"
%                    S(ng).pm              inverse grid X-spacing
%                    S(ng).pm_sum          interior accumulated "pm"
%                    S(ng).pm_diff         coarse minus averaged "pm"
%                    S(ng).pn              inverse grid Y-spacing
%                    S(ng).pn_sum          interior accumulated "pn"
%                    S(ng).pn_diff         coarse minus averaged "pn"
%                    S(ng).dndx            curvilinear metric d(n)/d(xi)
%                    S(ng).dndx_avg        interior averaged "dndx"
%                    S(ng).dndx_diff       coarse minus averaged "dndx"
%                    S(ng).dmde            curvilinear metric d(m)/d(eta)
%                    S(ng).dmde_avg        interior averaged "dmde"
%                    S(ng).dmde_diff       coarse minus averaged "dmde"
%                    S(ng).area            cell area (m^2)
%                    S(ng).area_sum        interior accumulated "area"
%                    S(ng).area_diff       coarse minus averaged "area"
%                    S(ng).volume          cell colunm volume (1E6 m^3)
%                    S(ng).volume_sum      interior accumulated "volume"
%                    S(ng).volume_diff     coarse minus averaged "volume"
%
%  The average values are similar to the fine-to-coarse operation in ROMS.
%  The interior finer grid values are averages at the location of the
%  coarser donor grid with coordinates (S(ng).xavg, S(ng).yavg).
%
%  The 2D linear indices, S(ng).dg_index, can be used to replace the
%  coarse donor grid values with the refined grid averaged values.  For
%  example:
%             index = S(rg).dg_index;     S(dg).f(index) = S(rg).f_avg
%
%  The difference fields correspond to the coarse minus the averaged finer
%  grid values.  Since we are using conservative quadratic interpolation,
%  the difference should be small.
%

% svn $Id$
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
% Set nested grid structure, G.
%--------------------------------------------------------------------------

if (~isstruct(Ginp))
  G = grids_structure(Ginp);
else
  G = Ginp;
end

Ngrids = length(G);

% Set ROMS grid size indices.

L = [G.Lm]+1;  Lp = L+1;
M = [G.Mm]+1;  Mp = M+1;

%--------------------------------------------------------------------------
% Set nested grids connectivity structure, C.
%--------------------------------------------------------------------------

if (~isstruct(Cinp))
  C = read_contact(Cinp);
else
  C = Cinp;
end

Ncontact = C.Ncontact;

% Set refined grid extraction coordinates from coarser donor grid.

Imin = double([C.grid.parent_Imin]);
Imax = double([C.grid.parent_Imax]);
Jmin = double([C.grid.parent_Jmin]);
Jmax = double([C.grid.parent_Jmax]);

%--------------------------------------------------------------------------
%  Average several refine grid fields to interior coarse grid cell:
%  fine-to-coarse step to check reversibility property.
%--------------------------------------------------------------------------

S(1:Ngrids) = struct('grid_name'    , [], 'spherical'    , [],          ...
                     'refine_factor', [],                               ...
                     'parent_Imin'  , [], 'parent_Imax'  , [],          ...
                     'parent_Jmin'  , [], 'parent_Jmax'  , [],          ...
                     'x'            , [], 'y'            , [],          ...
                     'Idg_avg'      , [], 'Jdg_avg'      , [],          ...
                     'dg_index'     , [],                               ...
                     'x_avg'        , [], 'y_avg'        , [],          ...
                     'x_perimeter'  , [], 'y_perimeter'  , [],          ...
                     'f'            , [],                               ...
                     'f_avg'        , [], 'f_diff'       , [],          ...
                     'h'            , [],                               ...
                     'h_avg'        , [], 'h_diff'       , [],          ...
                     'angle'        , [],                               ...
                     'angle_avg'    , [], 'angle_diff'   , [],          ...
                     'pm'           , [],                               ...
                     'pm_sum'       , [], 'pm_diff'      , [],          ...
                     'pn'           , [],                               ...
                     'pn_sum'       , [], 'pn_diff'      , [],          ...
                     'dndx'         , [],                               ...
                     'dndx_avg'     , [], 'dndx_diff'    , [],          ...
                     'dmde'         , [],                               ...
                     'dmde_avg'     , [], 'dmde_diff'    , [],          ...
                     'area'         , [],                               ...
                     'area_sum'     , [], 'area_diff'    , [],          ...
                     'volume'       , [],                               ...
                     'volume_sum'   , [], 'volume_diff'  , []);

AreaAvg = zeros([1 Ngrids]);

for ng = 1:Ngrids

  S(ng).grid_name     = G(ng).grid_name;
  S(ng).spherical     = G(ng).spherical;
  S(ng).refine_factor = G(ng).refine_factor;
  S(ng).parent_Imin   = Imin(ng); 
  S(ng).parent_Imax   = Imax(ng); 
  S(ng).parent_Jmin   = Jmin(ng); 
  S(ng).parent_Jmax   = Jmax(ng); 

  IstrP = 1;    IendP = G(ng).Lm+1;
  JstrP = 1;    JendP = G(ng).Mm+1;

  if (S(ng).spherical)
    X    = G(ng).lon_rho;
    Y    = G(ng).lat_rho;
    Xbox = [squeeze(G(ng).lon_psi(IstrP:IendP,JstrP));                  ...
            squeeze(G(ng).lon_psi(IendP,JstrP+1:JendP))';               ...
            squeeze(flipud(G(ng).lon_psi(IstrP:IendP-1,JendP)));        ...
            squeeze(fliplr(G(ng).lon_psi(IstrP,JstrP:JendP-1)))'];
    Ybox = [squeeze(G(ng).lat_psi(IstrP:IendP,JstrP));                  ...
            squeeze(G(ng).lat_psi(IendP,JstrP+1:JendP))';               ...
            squeeze(flipud(G(ng).lat_psi(IstrP:IendP-1,JendP)));        ...
            squeeze(fliplr(G(ng).lat_psi(IstrP,JstrP:JendP-1)))'];  
  else
    X    = G(ng).x_rho;
    Y    = G(ng).y_rho;
    Xbox = [squeeze(G(ng).x_psi(IstrP:IendP,JstrP));                    ...
            squeeze(G(ng).x_psi(IendP,JstrP+1:JendP))';                 ...
            squeeze(flipud(G(ng).x_psi(IstrP:IendP-1,JendP)));          ...
            squeeze(fliplr(G(ng).x_psi(IstrP,JstrP:JendP-1)))'];
    Ybox = [squeeze(G(ng).y_psi(IstrP:IendP,JstrP));                    ...
            squeeze(G(ng).y_psi(IendP,JstrP+1:JendP))';                 ...
            squeeze(flipud(G(ng).y_psi(IstrP:IendP-1,JendP)));          ...
            squeeze(fliplr(G(ng).y_psi(IstrP,JstrP:JendP-1)))'];
  end    

  S(ng).x           = X;
  S(ng).y           = Y;
  S(ng).x_perimeter = Xbox;
  S(ng).y_perimeter = Ybox;

  S(ng).f           = G(ng).f;
  S(ng).h           = G(ng).h;
  S(ng).angle       = G(ng).angle;
  S(ng).pm          = G(ng).pm;
  S(ng).pn          = G(ng).pn;
  S(ng).area        = (1./S(ng).pm) .* (1./S(ng).pn);
  S(ng).volume      = S(ng).area .* S(ng).h;
  
  if (~isempty(G(ng).dndx))
    S(ng).dndx      = G(ng).dndx;  
  else
    S(ng).dndx      = zeros(size(G(ng).h));  
  end

  if (~isempty(G(ng).dmde))
    S(ng).dmde      = G(ng).dmde;  
  else
    S(ng).dmde      = zeros(size(G(ng).h));  
  end
  
  rfactor = S(ng).refine_factor;

  if (rfactor > 0)
    [S(ng).f_avg, Xavg, Yavg] = refine_avg(X, Y, S(ng).f, rfactor);
    S(ng).x_avg = Xavg;
    S(ng).y_avg = Yavg;

    S(ng).h_avg      = refine_avg(X, Y, S(ng).h,      rfactor);

    S(ng).angle_avg  = refine_avg(X, Y, S(ng).angle,  rfactor);

    S(ng).dndx_avg   = refine_avg(X, Y, S(ng).dndx,   rfactor);
   
    S(ng).dmde_avg   = refine_avg(X, Y, S(ng).dmde,   rfactor);

    S(ng).pm_sum     = refine_sum(X, Y, S(ng).pm,     rfactor, 'pm');

    S(ng).pn_sum     = refine_sum(X, Y, S(ng).pn,     rfactor, 'pn');

    S(ng).area_sum   = refine_sum(X, Y, S(ng).area,   rfactor, 'area');

    S(ng).volume_sum = refine_sum(X, Y, S(ng).volume, rfactor, 'volume');
  end  
  
  AreaAvg(ng) = mean(S(ng).area(:));

end

% Compute difference fields between averaged and accumulated fields
% with coarser donor grid.

for cr=1:Ncontact
  dg = double(C.contact(cr).donor_grid);
  rg = double(C.contact(cr).receiver_grid);
  rfactor = S(rg).refine_factor;
  if ((rfactor > 0) && AreaAvg(dg) > AreaAvg(rg))
    Imid = 2+(rfactor-1)/2:rfactor:Lp(rg);
    Jmid = 2+(rfactor-1)/2:rfactor:Mp(rg);

    delta = 1.0/rfactor;
    half  = 0.5*delta;

    XpF = (Imin(rg):delta:Imax(rg));                 % PSI-points
    YpF = (Jmin(rg):delta:Jmax(rg));                 % PSI-points
    XrF = [XpF(1)-half XpF+half];                    % RHO-points
    YrF = [YpF(1)-half YpF+half];                    % RHO-points  

    Idg = fix(XrF(Imid) - 0.5) + 1;
    Jdg = fix(YrF(Jmid) - 0.5) + 1;
    
    [Yrg,Xrg] = meshgrid(Jdg,Idg);
    dg_index  = sub2ind([Lp(dg), Mp(dg)], Xrg, Yrg);
    
    S(rg).Idg_avg     = Idg;
    S(rg).Jdg_avg     = Jdg;
    S(rg).dg_index    = dg_index;
    
    S(rg).f_diff      = S(dg).f     (dg_index) - S(rg).f_avg;
    S(rg).h_diff      = S(dg).h     (dg_index) - S(rg).h_avg;
    S(rg).angle_diff  = S(dg).angle (dg_index) - S(rg).angle_avg;
    S(rg).dndx_diff   = S(dg).dndx  (dg_index) - S(rg).dndx_avg;
    S(rg).dmde_diff   = S(dg).dmde  (dg_index) - S(rg).dmde_avg;

    S(rg).pm_diff     = S(dg).pm    (dg_index) - S(rg).pm_sum;
    S(rg).pn_diff     = S(dg).pn    (dg_index) - S(rg).pn_sum;
    S(rg).area_diff   = S(dg).area  (dg_index) - S(rg).area_sum;
    S(rg).volume_diff = S(dg).volume(dg_index) - S(rg).volume_sum;
  end
end

%--------------------------------------------------------------------------
% Plot various fields.
%--------------------------------------------------------------------------

if (Lplot)
  cst = 'm-';
  box = 'r-';

  Mversion = version('-release');
  Vyear    = sscanf(Mversion, '%i');

  figure;

  if (strcmp(Mversion,'2014b') || Vyear > 2014)
    h = gcf;
    nf = h.Number;
  else
    ag = findobj;                                   % all graphical objects
    nf = max(ag(ag == fix(ag)));                    % number of figures
  end

  h_diff      = cell(1, Ngrids-1);
  f_diff      = cell(1, Ngrids-1);
  angle_diff  = cell(1, Ngrids-1);
  pm_diff     = cell(1, Ngrids-1);
  pn_diff     = cell(1, Ngrids-1);
  dndx_diff   = cell(1, Ngrids-1);
  dmde_diff   = cell(1, Ngrids-1);
  area_diff   = cell(1, Ngrids-1);
  volume_diff = cell(1, Ngrids-1);

  for ng=1:Ngrids
    ic = nf-1;

    if (S(ng).spherical)
      Xavg = S(ng).x_avg;
      Yavg = S(ng).y_avg;
      Xperimeter = S(ng).x_perimeter;
      Yperimeter = S(ng).y_perimeter;
    else
      Xavg = 0.001 .* S(ng).x_avg;                  % km
      Yavg = 0.001 .* S(ng).y_avg;
      Xperimeter = 0.001 .* S(ng).x_perimeter;
      Yperimeter = 0.001 .* S(ng).y_perimeter;
    end

    if (ng == 1)    
      got.coast = false;
      if (S(ng).spherical)
        Xmin = min(S(ng).x(:));
        Xmax = max(S(ng).x(:));
        Ymin = min(S(ng).y(:));
        Ymax = max(S(ng).y(:));

        if (isfield(G(ng),'lon_coast') && isfield(G(ng),'lat_coast'))
          Clon = G(ng).lon_coast;
          Clat = G(ng).lat_coast;
          got.coast = true;
        end
      else
        Xmin = min(0.001 .* S(ng).x(:));
        Xmax = max(0.001 .* S(ng).x(:));
        Ymin = min(0.001 .* S(ng).y(:));
        Ymax = max(0.001 .* S(ng).y(:));
      end

      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Bathymetry Differences: Fine-to-Coarse');

      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Coriolis Differences: Fine-to-Coarse');
      
      if (G(ng).curvilinear)
        ic = ic + 1;
        figure(ic);
        plot(Xperimeter, Yperimeter, box);
        axis([Xmin Xmax Ymin Ymax]); grid on;
        hold on;
        title('Curvilinear Angle Differences: Fine-to-Coarse');
      end
 
      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Inverse Grid X-Spacing Differences: Fine-to-Coarse');

      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Inverse Grid Y-Spacing Differences: Fine-to-Coarse');
      
      if (G(ng).curvilinear)
        ic = ic + 1;
        figure(ic);
        plot(Xperimeter, Yperimeter, box);
        axis([Xmin Xmax Ymin Ymax]); grid on;
        hold on;
        title('Curvilinear Metric DNDX Differences: Fine-to-Coarse');

        ic = ic + 1;
        figure(ic);
        plot(Xperimeter, Yperimeter, box);
        axis([Xmin Xmax Ymin Ymax]); grid on;
        hold on;
        title('Curvilinear Metric DMDE Differences: Fine-to-Coarse');
      end

      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Area (squared meters) Differences: Fine-to-Coarse');

      ic = ic + 1;
      figure(ic);
      plot(Xperimeter, Yperimeter, box);
      axis([Xmin Xmax Ymin Ymax]); grid on;
      hold on;
      title('Volume (millions cubic meters) Differences: Fine-to-Coarse');
    
    else

      ic = ic + 1;
      figure(ic);      
      pcolorjw(Xavg, Yavg, S(ng).h_diff); shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      Vmin = min(S(ng).h_diff(:));
      Vmax = max(S(ng).h_diff(:));
      string = ['Grid = ', num2str(ng), ',', blanks(4),                 ...
                'Min = ', num2str(Vmin, '% 0.5e'), blanks(4),           ...
                'Max = ', num2str(Vmax, '% 0.5e')];
      h_diff(ng-1) = {string};
      if (ng == Ngrids)
        xlabel(h_diff);
        if (got.coast)
          plot(Clon, Clat, cst);
        end
      end

      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xavg, Yavg, S(ng).f_diff); shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      Vmin = min(S(ng).f_diff(:));
      Vmax = max(S(ng).f_diff(:));
      string = ['Grid = ', num2str(ng), ',', blanks(4),                 ...
                'Min = ', num2str(Vmin, '% 0.5e'), blanks(4),           ...
                'Max = ', num2str(Vmax, '% 0.5e')];
      f_diff(ng-1) = {string};
      if (ng == Ngrids)
        xlabel(f_diff);
        if (got.coast)
          plot(Clon, Clat, cst);
        end
      end
      
      if (G(ng).curvilinear)
        ic = ic + 1;
        figure(ic);      
        pcolorjw (Xavg, Yavg, S(ng).angle_diff); shading flat; colorbar;
        plot(Xperimeter, Yperimeter, box);
        Vmin = min(S(ng).angle_diff(:));
        Vmax = max(S(ng).angle_diff(:));
        string = ['Grid = ', num2str(ng), ',', blanks(4),               ...
                  'Min = ', num2str(Vmin, '% 0.5e'), blanks(4),         ...
                  'Max = ', num2str(Vmax, '% 0.5e')];
        angle_diff(ng-1) = {string};
        if (ng == Ngrids)
          xlabel(angle_diff);
          if (got.coast)
            plot(Clon, Clat, cst);
          end
        end
      end

      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xavg, Yavg, S(ng).pm_diff); shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      Vmin = min(S(ng).pm_diff(:));
      Vmax = max(S(ng).pm_diff(:));
      string = ['Grid = ', num2str(ng), ',', blanks(4),                 ...
                'Min = ', num2str(Vmin, '% 0.5e'), blanks(4),           ...
                'Max = ', num2str(Vmax, '% 0.5e')];
      pm_diff(ng-1) = {string};
      if (ng == Ngrids)
        xlabel(pm_diff);
        if (got.coast)
          plot(Clon, Clat, cst);
        end
      end
    
      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xavg, Yavg, S(ng).pn_diff); shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      Vmin = min(S(ng).pn_diff(:));
      Vmax = max(S(ng).pn_diff(:));
      string = ['Grid = ', num2str(ng), ',', blanks(4),                 ...
                'Min = ', num2str(Vmin, '% 0.5e'), blanks(4),           ...
                'Max = ', num2str(Vmax, '% 0.5e')];
      pn_diff(ng-1) = {string};
      if (ng == Ngrids)
        xlabel(pn_diff);
        if (got.coast)
          plot(Clon, Clat, cst);
        end
      end
    
      if (G(ng).curvilinear)
        ic = ic + 1;
        figure(ic);      
        pcolorjw (Xavg, Yavg, S(ng).dndx_diff); shading flat; colorbar;
        plot(Xperimeter, Yperimeter, box);
        Vmin = min(S(ng).dndx_diff(:));
        Vmax = max(S(ng).dndx_diff(:));
        string = ['Grid = ', num2str(ng), ',', blanks(4),               ...
                  'Min = ', num2str(Vmin, '% 0.5e'), blanks(4),         ...
                  'Max = ', num2str(Vmax, '% 0.5e')];
        dndx_diff(ng-1) = {string};
        if (ng == Ngrids)
          xlabel(dndx_diff);
          if (got.coast)
            plot(Clon, Clat, cst);
          end
        end

        ic = ic + 1;
        figure(ic);      
        pcolorjw (Xavg, Yavg, S(ng).dmde_diff); shading flat; colorbar;
        plot(Xperimeter, Yperimeter, box);
        Vmin = min(S(ng).dmde_diff(:));
        Vmax = max(S(ng).dmde_diff(:));
        string = ['Grid = ', num2str(ng), ',', blanks(4),               ...
                  'Min = ', num2str(Vmin, '% 0.5e'), blanks(4),         ...
                  'Max = ', num2str(Vmax, '% 0.5e')];
        dmde_diff(ng-1) = {string};
        if (ng == Ngrids)
          xlabel(dmde_diff);
          if (got.coast)
            plot(Clon, Clat, cst);
          end
        end
      end

      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xavg, Yavg, S(ng).area_diff);
      shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      Vmin = min(S(ng).area_diff(:));
      Vmax = max(S(ng).area_diff(:));
      string = ['Grid = ', num2str(ng), ',', blanks(4),                 ...
                'Min = ', num2str(Vmin, '% 0.5e'), blanks(4),           ...
                'Max = ', num2str(Vmax, '% 0.5e')];
      area_diff(ng-1) = {string};
      if (ng == Ngrids)
        xlabel(area_diff);
        if (got.coast)
          plot(Clon, Clat, cst);
        end
      end
      
      ic = ic + 1;
      figure(ic);      
      pcolorjw (Xavg, Yavg, 1.0d-6 .* S(ng).volume_diff);
      shading flat; colorbar;
      plot(Xperimeter, Yperimeter, box);
      Vmin = min(1.0d-6 .* S(ng).volume_diff(:));
      Vmax = max(1.0d-6 .* S(ng).volume_diff(:));
      string = ['Grid = ', num2str(ng), ',', blanks(4),                 ...
                'Min = ', num2str(Vmin, '% 0.5e'), blanks(4),           ...
                'Max = ', num2str(Vmax, '% 0.5e')];
      volume_diff(ng-1) = {string};
      if (ng == Ngrids)
        xlabel(volume_diff);
        if (got.coast)
          plot(Clon, Clat, cst);
        end
      end
    end
  end
end  
      
return

%--------------------------------------------------------------------------

function [Favg, Xavg, Yavg] = refine_avg(X, Y, F, rfactor)

% This function averages requested field to donor coarse grid cells. This
% averaging is what ROMS does in subroutine "fine2coarse". The averages
% are computed at the finer grid interior points.
%
% On Input:
%
%    X          Finer grid X-coordinates (Cartesian or spherical; array)
%    Y          Finer grid Y-coordinates (Cartesian or spherical; array)
%    F          Finer field to process (array).
%    rfactor    Grid refinement factor  
%
% On Output:
%
%     Favg      Averaged field values (array)
%     Xavg      Averaged field X-coordinate (OPTIONAL; array)
%     Yavg      Averaged field Y-coordinate (OPTIONAL; array)
%
  
% Set averaged field coordinates.

[Im,Jm] = size(F);
  
half = (rfactor-1)/2;

Imid = 2+half:rfactor:Im;
Jmid = 2+half:rfactor:Jm;

if (nargout > 1)
  Xavg = X(Imid, Jmid);
  Yavg = Y(Imid, Jmid);
end

% Set averaged field.

Istr = 2:rfactor:Im-1;
Iend = Istr(2)-1:rfactor:Im;

Jstr = 2:rfactor:Jm-1;
Jend = Jstr(2)-1:rfactor:Jm;

Favg = nan([length(Imid) length(Jmid)]);
for j = 1:length(Jstr)
  for i = 1:length(Istr)
    values = F(Istr(i):Iend(i), Jstr(j):Jend(j));
    Favg(i,j) = mean(values(:));
  end
end

return

%--------------------------------------------------------------------------

function [Fsum, Xsum, Ysum] = refine_sum(X, Y, F, rfactor, vname)

% This function averages requested field to donor coarse grid cells. This
% averaging is what ROMS does in subroutine "fine2coarse". The averages
% are computed at the finer grid interior points.
%
% On Input:
%
%    X          Finer grid X-coordinates (Cartesian or spherical; array)
%    Y          Finer grid Y-coordinates (Cartesian or spherical; array)
%    F          Finer field to process (array).
%    rfactor    Grid refinement factor (3,5,7,9,...)  
%    vname      variable name (string)
%
% On Output:
%
%     Fsum      Accumulated field values (array)
%     Xsum      Accumulated field X-coordinate (OPTIONAL; array)
%     Ysum      Accumulated field Y-coordinate (OPTIONAL; array)
%
  
% Set accumulated field coordinates.

[Im,Jm] = size(F);

half = (rfactor-1)/2;

Imid = 2+half:rfactor:Im;
Jmid = 2+half:rfactor:Jm;

if (nargout > 1)
  Xsum = X(Imid, Jmid);
  Ysum = Y(Imid, Jmid);
end

% Set accumulated field.

Istr = 2:rfactor:Im-1;
Iend = Istr(2)-1:rfactor:Im;

Jstr = 2:rfactor:Jm-1;
Jend = Jstr(2)-1:rfactor:Jm;

switch (vname)
  case {'pm'}
    dx = nan([length(Imid) length(Jmid)]);
    for j = 1:length(Jstr)
      Javg = Jmid(j)-half:1:Jmid(j)+half;
      Nvalues = length(Javg);
      values  = zeros([1 Nvalues]);
      for i = 1:length(Istr)
        for jc = 1:Nvalues
          values(jc) = sum(1 ./ squeeze(F(Istr(i):Iend(i), Javg(jc))));
        end
        dx(i,j) = mean(values(:));
      end
    end
    Fsum = 1 ./ dx;
  case {'pn'} 
    dy = nan([length(Imid) length(Jmid)]);
    for i = 1:length(Istr)
      Iavg = Imid(i)-half:1:Imid(i)+half;
      Nvalues = length(Iavg);
      values  = zeros([1 Nvalues]);
      for j = 1:length(Jstr)
        for ic = 1:Nvalues
          values(ic) = sum(1 ./ squeeze(F(Iavg(ic), Jstr(j):Jend(j))));
        end
        dy(i,j) = mean(values(:));
      end
    end 
    Fsum = 1 ./ dy;
 otherwise
    Fsum = nan([length(Imid) length(Jmid)]);
    for j = 1:length(Jstr)
      for i = 1:length(Istr)
        values = F(Istr(i):Iend(i), Jstr(j):Jend(j));
        Fsum(i,j) = sum(values(:));
      end
    end
end  

return
