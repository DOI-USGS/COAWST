function S = grid_perimeter(G)

%
% GRID_PERIMETER:  Sets Nested Grids Perimeters and Boundary Edges
%
% S = grid_perimeter(G)
%
% This function creates a structure containing information about nested
% grids perimeters, boundary edges, and other parameters.
%
% On Input:
%
%    G          Information grids structure (1 x Ngrid struct array)
%
%                 G(ng) = get_roms_grid ( char(Gnames(ng)) )
%
% On Output:
%
%    S          Nested grids information structure (struct array)
%
%
% This function adds the following fields to the output structure:
%
%    S.Ngrids                         - Number of nested grids
%    S.Ncontact                       - Number of contact regions
%    S.nLweights = 4                  - Number of linear weights
%    S.nQweights = 9                  - Number of quadratic weights
%    S.Ndatum = 0                     - Total number of contact points
%
%    S.western_edge  = 1              - Western  boundary edge index
%    S.southern_edge = 2              - Southern boundary edge index
%    S.eastern_edge  = 3              - Eastern  boundary edge index
%    S.northern_edge = 4              - Northern boundary edge index
%
%    S.spherical                      - Spherical switch
%
%    S.grid(ng).filename              - Grid NetCDF file name
%
%    S.grid(ng).Lp                    - Number of I-points (RHO)
%    S.grid(ng).Mp                    - Number of J-points (RHO)
%    S.grid(ng).L                     - Number of I-points (PSI)
%    S.grid(ng).M                     - Number of J-points (PSI)
%
%    S.grid(ng).refine_factor         - Refinement factor (0,3,5,7)
%
%    S.grid(ng).XI_psi (:,:)          - ROMS XI-coordinates  (PSI)
%    S.grid(ng).ETA_psi(:,:)          - ROMS ETA-coordinates (PSI)
%    S.grid(ng).XI_rho (:,:)          - ROMS XI-coordinates  (RHO)
%    S.grid(ng).ETA_rho(:,:)          - ROMS ETA-coordinates (RHO)
%    S.grid(ng).XI_u   (:,:)          - ROMS XI-coordinates  (U)
%    S.grid(ng).ETA_u  (:,:)          - ROMS ETA-coordinates (U)
%    S.grid(ng).XI_v   (:,:)          - ROMS XI-coordinates  (V)
%    S.grid(ng).ETA_v  (:,:)          - ROMS ETA-coordinates (V)
%
%    S.grid(ng).I_psi(:,:)            - ROMS I-indices at PSI-points
%    S.grid(ng).J_psi(:,:)            - ROMS J-indices at PSI-points
%    S.grid(ng).I_rho(:,:)            - ROMS I-indices at RHO-points
%    S.grid(ng).J_rho(:,:)            - ROMS J-indices at RHO-points
%    S.grid(ng).I_u  (:,:)            - ROMS I-indices at U-points
%    S.grid(ng).J_u  (:,:)            - ROMS J-indices at U-points
%    S.grid(ng).I_v  (:,:)            - ROMS I-indices at V-points
%    S.grid(ng).J_v  (:,:)            - ROMS I-indices at V-points
%
%    S.grid(ng).perimeter.X_psi(:)    - Perimeter X-coordinates (PSI)
%    S.grid(ng).perimeter.Y_psi(:)    - Perimeter Y-coordinates (PSI)
%    S.grid(ng).perimeter.X_rho(:)    - Perimeter X-coordinates (RHO)
%    S.grid(ng).perimeter.Y_rho(:)    - Perimeter Y-coordinates (RHO)
%    S.grid(ng).perimeter.X_u(:)      - Perimeter X-coordinates (U)
%    S.grid(ng).perimeter.Y_u(:)      - Perimeter Y-coordinates (U)
%    S.grid(ng).perimeter.X_v(:)      - Perimeter X-coordinates (V)
%    S.grid(ng).perimeter.Y_v(:)      - Perimeter Y-coordinates (V)
%    S.grid(ng).perimeter.X_uv(:)     - Perimeter X-coordinates (U-V)
%    S.grid(ng).perimeter.Y_uv(:)     - Perimeter Y-coordinates (U-V)
%                                       (counterclockwise)
%
%    S.grid(ng).corners.index(:)      - Corners linear IJ-index
%    S.grid(ng).corners.X(:)          - Corners X-coordinates (PSI)
%    S.grid(ng).corners.Y(:)          - Corners Y-coordinates (PSI)
%    S.grid(ng).corners.I(:)          - Corners I-indices (PSI)
%    S.grid(ng).corners.J(:)          - Corners J-indices (PSI)
%
%    S.grid(ng).boundary(ib).index(:) - Boundary linear IJ-index
%    S.grid(ng).boundary(ib).X(:)     - Boundary X-coordinates (PSI)
%    S.grid(ng).boundary(ib).Y(:)     - Boundary Y-coordinates (PSI)
%    S.grid(ng).boundary(ib).Xuv(:)   - Boundary X-coordinates (U,V)
%    S.grid(ng).boundary(ib).Yuv(:)   - Boundary Y-coordinates (U,V)
%

% svn $Id: grid_perimeter.m 738 2014-10-14 21:49:14Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

S.Ngrids    = length(G);
S.Ncontact  = (S.Ngrids-1)*2;
S.nLweights = 4;
S.nQweights = 9;
S.Ndatum    = 0;

S.western_edge  = 1;
S.southern_edge = 2;
S.eastern_edge  = 3;
S.northern_edge = 4;

S.spherical = G(1).spherical;

for ng=1:S.Ngrids,
  S.grid(ng).filename = G(ng).grid_name;
  
  S.grid(ng).Lp = G(ng).Lm+2;
  S.grid(ng).Mp = G(ng).Mm+2;

  S.grid(ng).L  = G(ng).Lm+1;
  S.grid(ng).M  = G(ng).Mm+1;
end

% Get grid refinement factor, if any. Otherwise, set to zero.

for ng=1:S.Ngrids,
  S.grid(ng).refine_factor = 0;
  if (isfield(G(ng),'refine_factor')),
    if (~isempty(G(ng).refine_factor)),
      S.grid(ng).refine_factor = G(ng).refine_factor;
    end
  end
end

%--------------------------------------------------------------------------
% Set (XI,ETA) fractional coordinates.
%--------------------------------------------------------------------------

for ng=1:S.Ngrids,
  Lp = S.grid(ng).Lp;   Mp = S.grid(ng).Mp;
  L  = S.grid(ng).L;    M  = S.grid(ng).M;

  [Yp, Xp] = meshgrid(1.0:1:M     , 1.0:1:L     );      % PSI-points
  S.grid(ng).XI_psi  = Xp;
  S.grid(ng).ETA_psi = Yp;

  [Yr, Xr] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);      % RHO-points
  S.grid(ng).XI_rho  = Xr;
  S.grid(ng).ETA_rho = Yr;
  
  [Yu, Xu] = meshgrid(0.5:1:Mp-0.5, 1.0:1:L     );      % U-points
  S.grid(ng).XI_u    = Xu;
  S.grid(ng).ETA_u   = Yu;

  [Yv, Xv] = meshgrid(1.0:1:M     , 0.5:1:Lp-0.5);      % V-points
  S.grid(ng).XI_v    = Xv;
  S.grid(ng).ETA_v   = Yv;

% Grid perimeters in (XI,ETA) coordinates (counterclockwise from south).

  IstrP = 1;      IstrR = 1;        IstrU = 1;        IstrV = 1;
  IendP = L;      IendR = L+1;      IendU = L;        IendV = L+1;
  JstrP = 1;      JstrR = 1;        JstrU = 1;        JstrV = 1;
  JendP = M;      JendR = M+1;      JendU = M+1;      JendV = M;

% PSI-points.

  Xbox = [squeeze(Xp(IstrP:IendP,JstrP));                               ...
          squeeze(Xp(IendP,JstrP+1:JendP))';                            ...
          squeeze(flipud(Xp(IstrP:IendP-1,JendP)));                     ...
          squeeze(fliplr(Xp(IstrP,JstrP:JendP-1)))'];
    
  Ybox = [squeeze(Yp(IstrP:IendP,JstrP));                               ...
          squeeze(Yp(IendP,JstrP+1:JendP))';                            ...
          squeeze(flipud(Yp(IstrP:IendP-1,JendP)));                     ...
          squeeze(fliplr(Yp(IstrP,JstrP:JendP-1)))'];

  S.grid(ng).perimeter.XI_psi  = Xbox;
  S.grid(ng).perimeter.ETA_psi = Ybox;

% RHO-points.

  Xbox = [squeeze(Xr(IstrR:IendR,JstrR));                               ...
          squeeze(Xr(IendR,JstrR+1:JendR))';                            ...
          squeeze(flipud(Xr(IstrR:IendR-1,JendR)));                     ...
          squeeze(fliplr(Xr(IstrR,JstrR:JendR-1)))'];
    
  Ybox = [squeeze(Yr(IstrR:IendR,JstrR));                               ...
          squeeze(Yr(IendR,JstrR+1:JendR))';                            ...
          squeeze(flipud(Yr(IstrR:IendR-1,JendR)));                     ...
          squeeze(fliplr(Yr(IstrR,JstrR:JendR-1)))'];

  S.grid(ng).perimeter.XI_rho  = Xbox;
  S.grid(ng).perimeter.ETA_rho = Ybox;

% U-points.

  Xbox = [squeeze(Xu(IstrU:IendU,JstrU));                               ...
          squeeze(Xu(IendU,JstrU+1:JendU))';                            ...
          squeeze(flipud(Xu(IstrU:IendU-1,JendU)));                     ...
          squeeze(fliplr(Xu(IstrU,JstrU:JendU-1)))'];
    
  Ybox = [squeeze(Yu(IstrU:IendU,JstrU));                               ...
          squeeze(Yu(IendU,JstrU+1:JendU))';                            ...
          squeeze(flipud(Yu(IstrU:IendU-1,JendU)));                     ...
          squeeze(fliplr(Yu(IstrU,JstrU:JendU-1)))'];

  S.grid(ng).perimeter.XI_u  = Xbox;
  S.grid(ng).perimeter.ETA_u = Ybox;

% V-points.

  Xbox = [squeeze(Xv(IstrV:IendV,JstrV));                               ...
          squeeze(Xv(IendV,JstrV+1:JendV))';                            ...
          squeeze(flipud(Xv(IstrV:IendV-1,JendV)));                     ...
          squeeze(fliplr(Xv(IstrV,JstrV:JendV-1)))'];
    
  Ybox = [squeeze(Yv(IstrV:IendV,JstrV));                               ...
          squeeze(Yv(IendV,JstrV+1:JendV))';                            ...
          squeeze(flipud(Yv(IstrV:IendV-1,JendV)));                     ...
          squeeze(fliplr(Yv(IstrV,JstrV:JendV-1)))'];

    S.grid(ng).perimeter.XI_v  = Xbox;
    S.grid(ng).perimeter.ETA_v = Ybox;

% UV-points.

  Xbox = [squeeze(Xv(IstrV+1:IendV-1,JstrV));                           ...
          squeeze(Xu(IendU,JstrU+1:JendU-1))';                          ...
          squeeze(flipud(Xv(IstrV+1:IendV-1,JendV)));                   ...
          squeeze(fliplr(Xu(IstrU,JstrU+1:JendU-1)))'];
    
  Ybox = [squeeze(Yv(IstrV+1:IendV-1,JstrV));                           ...
          squeeze(Yu(IendU,JstrU+1:JendU-1))';                          ...
          squeeze(flipud(Yv(IstrV+1:IendV-1,JendV)));                   ...
          squeeze(fliplr(Yu(IstrU,JstrU+1:JendU-1)))'];

  S.grid(ng).perimeter.XI_uv  = Xbox;
  S.grid(ng).perimeter.ETA_uv = Ybox;

end

%--------------------------------------------------------------------------
% Spherical grids: set grid indices, perimeters, corners, and boundary
%                  edges.
%--------------------------------------------------------------------------

if (S.spherical),

  for ng=1:S.Ngrids,
    Im = S.grid(ng).L;
    Jm = S.grid(ng).M;
    
    IstrP = 1;      IstrR = 1;        IstrU = 1;        IstrV = 1;
    IendP = Im;     IendR = Im+1;     IendU = Im;       IendV = Im+1;
    JstrP = 1;      JstrR = 1;        JstrU = 1;        JstrV = 1;
    JendP = Jm;     JendR = Jm+1;     JendU = Jm+1;     JendV = Jm;
    
% C-type variables ROMS indices. 

    [S.grid(ng).J_psi, S.grid(ng).I_psi] = meshgrid(1:Jm, 1:Im);
    [S.grid(ng).J_rho, S.grid(ng).I_rho] = meshgrid(0:Jm, 0:Im);
    [S.grid(ng).J_u  , S.grid(ng).I_u  ] = meshgrid(0:Jm, 1:Im);
    [S.grid(ng).J_v  , S.grid(ng).I_v  ] = meshgrid(1:Jm, 0:Im);

% Grid perimeter at PSI-points (counterclockwise from south). This is the
% physical grid perimeter.

    Xbox = [squeeze(G(ng).lon_psi(IstrP:IendP,JstrP));                  ...
            squeeze(G(ng).lon_psi(IendP,JstrP+1:JendP))';               ...
            squeeze(flipud(G(ng).lon_psi(IstrP:IendP-1,JendP)));        ...
            squeeze(fliplr(G(ng).lon_psi(IstrP,JstrP:JendP-1)))'];
    
    Ybox = [squeeze(G(ng).lat_psi(IstrP:IendP,JstrP));                  ...
            squeeze(G(ng).lat_psi(IendP,JstrP+1:JendP))';               ...
            squeeze(flipud(G(ng).lat_psi(IstrP:IendP-1,JendP)));        ...
            squeeze(fliplr(G(ng).lat_psi(IstrP,JstrP:JendP-1)))'];

    S.grid(ng).perimeter.X_psi = Xbox;
    S.grid(ng).perimeter.Y_psi = Ybox;

% Grid perimeter at RHO-points (counterclockwise from south). Needed for
% lateral boundary condition switch at RHO-points.

    Xbox = [squeeze(G(ng).lon_rho(IstrR:IendR,JstrR));                  ...
            squeeze(G(ng).lon_rho(IendR,JstrR+1:JendR))';               ...
            squeeze(flipud(G(ng).lon_rho(IstrR:IendR-1,JendR)));        ...
            squeeze(fliplr(G(ng).lon_rho(IstrR,JstrR:JendR-1)))'];
    
    Ybox = [squeeze(G(ng).lat_rho(IstrR:IendR,JstrR));                  ...
            squeeze(G(ng).lat_rho(IendR,JstrR+1:JendR))';               ...
            squeeze(flipud(G(ng).lat_rho(IstrR:IendR-1,JendR)));        ...
            squeeze(fliplr(G(ng).lat_rho(IstrR,JstrR:JendR-1)))'];

    S.grid(ng).perimeter.X_rho = Xbox;
    S.grid(ng).perimeter.Y_rho = Ybox;

% Grid perimeter at U-points (counterclockwise from south). Needed for
% lateral boundary condition switch at U-points.

    Xbox = [squeeze(G(ng).lon_u(IstrU:IendU,JstrU));                    ...
            squeeze(G(ng).lon_u(IendU,JstrU+1:JendU))';                 ...
            squeeze(flipud(G(ng).lon_u(IstrU:IendU-1,JendU)));          ...
            squeeze(fliplr(G(ng).lon_u(IstrU,JstrU:JendU-1)))'];
    
    Ybox = [squeeze(G(ng).lat_u(IstrU:IendU,JstrU));                    ...
            squeeze(G(ng).lat_u(IendU,JstrU+1:JendU))';                 ...
            squeeze(flipud(G(ng).lat_u(IstrU:IendU-1,JendU)));          ...
            squeeze(fliplr(G(ng).lat_u(IstrU,JstrU:JendU-1)))'];

    S.grid(ng).perimeter.X_u = Xbox;
    S.grid(ng).perimeter.Y_u = Ybox;

% Grid perimeter at V-points (counterclockwise from south). Needed for
% lateral boundary condition switch at V-points.

    Xbox = [squeeze(G(ng).lon_v(IstrV:IendV,JstrV));                    ...
            squeeze(G(ng).lon_v(IendV,JstrV+1:JendV))';                 ...
            squeeze(flipud(G(ng).lon_v(IstrV:IendV-1,JendV)));          ...
            squeeze(fliplr(G(ng).lon_v(IstrV,JstrV:JendV-1)))'];
    
    Ybox = [squeeze(G(ng).lat_v(IstrV:IendV,JstrV));                    ...
            squeeze(G(ng).lat_v(IendV,JstrV+1:JendV))';                 ...
            squeeze(flipud(G(ng).lat_v(IstrV:IendV-1,JendV)));          ...
            squeeze(fliplr(G(ng).lat_v(IstrV,JstrV:JendV-1)))'];

    S.grid(ng).perimeter.X_v = Xbox;
    S.grid(ng).perimeter.Y_v = Ybox;
    
% Grid perimeter at UV-points (counterclockwise from south). Needed for
% lateral boundary condition switch at U- and V-points. It is more
% accurate to use this perimeter with inpolygon.

    Xbox = [squeeze(G(ng).lon_v(IstrV+1:IendV-1,JstrV));                ...
            squeeze(G(ng).lon_u(IendU,JstrU+1:JendU-1))';               ...
            squeeze(flipud(G(ng).lon_v(IstrV+1:IendV-1,JendV)));        ...
            squeeze(fliplr(G(ng).lon_u(IstrU,JstrU+1:JendU-1)))'];
    
    Ybox = [squeeze(G(ng).lat_v(IstrV+1:IendV-1,JstrV));                ...
            squeeze(G(ng).lat_u(IendU,JstrU+1:JendU-1))';               ...
            squeeze(flipud(G(ng).lat_v(IstrV+1:IendV-1,JendV)));        ...
            squeeze(fliplr(G(ng).lat_u(IstrU,JstrU+1:JendU-1)))'];

    S.grid(ng).perimeter.X_uv = Xbox;
    S.grid(ng).perimeter.Y_uv = Ybox;
    
% Grid domain corners (PSI-points).

    Icorners = false([Im Jm]);
    Icorners(1 ,1 ) = true;
    Icorners(Im,1 ) = true;
    Icorners(Im,Jm) = true;
    Icorners(1 ,Jm) = true;
    Icorners = find(Icorners(:) == true);

    S.grid(ng).corners.index = Icorners;
    S.grid(ng).corners.X     = G(ng).lon_psi(Icorners);
    S.grid(ng).corners.Y     = G(ng).lat_psi(Icorners);
    S.grid(ng).corners.I     = S.grid(ng).I_psi(Icorners);
    S.grid(ng).corners.J     = S.grid(ng).J_psi(Icorners);

% Boundary edges at PSI-points, excluding corners. 

    Iwest  = false([Im Jm]);    Iwest (1 ,2:Jm-1) = true;
    Isouth = false([Im Jm]);    Isouth(2:Im-1, 1) = true;
    Ieast  = false([Im Jm]);    Ieast (Im,2:Jm-1) = true;
    Inorth = false([Im Jm]);    Inorth(2:Im-1,Jm) = true;

    B(S.western_edge ).index = find(Iwest (:) == true);
    B(S.southern_edge).index = find(Isouth(:) == true);
    B(S.eastern_edge ).index = find(Ieast (:) == true);
    B(S.northern_edge).index = find(Inorth(:) == true);

    for ib=1:4,
      S.grid(ng).boundary(ib).index = B(ib).index;
      S.grid(ng).boundary(ib).X     = G(ng).lon_psi(B(ib).index);
      S.grid(ng).boundary(ib).Y     = G(ng).lat_psi(B(ib).index);
    end
  
% Boundary edges at U-points and V-points: Momentum physical boundary. 

    Iwest  = false([Im Jm+1]);  Iwest (1 ,2:Jm) = true;
    Ieast  = false([Im Jm+1]);  Ieast (Im,2:Jm) = true;

    Isouth = false([Im+1 Jm]);  Isouth(2:Im, 1) = true;
    Inorth = false([Im+1 Jm]);  Inorth(2:Im,Jm) = true;

    B(1).index = find(Iwest (:) == true);
    B(2).index = find(Isouth(:) == true);
    B(3).index = find(Ieast (:) == true);
    B(4).index = find(Inorth(:) == true);

    S.grid(ng).boundary(1).Xuv = G(ng).lon_u(B(1).index);
    S.grid(ng).boundary(2).Xuv = G(ng).lon_v(B(2).index);
    S.grid(ng).boundary(3).Xuv = G(ng).lon_u(B(3).index);
    S.grid(ng).boundary(4).Xuv = G(ng).lon_v(B(4).index);

    S.grid(ng).boundary(1).Yuv = G(ng).lat_u(B(1).index);
    S.grid(ng).boundary(2).Yuv = G(ng).lat_v(B(2).index);
    S.grid(ng).boundary(3).Yuv = G(ng).lat_u(B(3).index);
    S.grid(ng).boundary(4).Yuv = G(ng).lat_v(B(4).index);
  
  end
end

%--------------------------------------------------------------------------
% Cartesian grids: set grid indices, perimeters, corners, and boundary
%                  edges.
%--------------------------------------------------------------------------

if (~S.spherical),

  for ng=1:S.Ngrids,
    Im = S.grid(ng).L;
    Jm = S.grid(ng).M;
    
    IstrP = 1;      IstrR = 1;        IstrU = 1;        IstrV = 1;
    IendP = Im;     IendR = Im+1;     IendU = Im;       IendV = Im+1;
    JstrP = 1;      JstrR = 1;        JstrU = 1;        JstrV = 1;
    JendP = Jm;     JendR = Jm+1;     JendU = Jm+1;     JendV = Jm;
    
% C-type variables ROMS indices. 

    [S.grid(ng).J_psi, S.grid(ng).I_psi] = meshgrid(1:Jm, 1:Im);
    [S.grid(ng).J_rho, S.grid(ng).I_rho] = meshgrid(0:Jm, 0:Im);
    [S.grid(ng).J_u  , S.grid(ng).I_u  ] = meshgrid(0:Jm, 1:Im);
    [S.grid(ng).J_v  , S.grid(ng).I_v  ] = meshgrid(1:Jm, 0:Im);

% Grid perimeter at PSI-points (counterclockwise from south). This is the
% physical grid perimeter.

    Xbox = [squeeze(G(ng).x_psi(IstrP:IendP,JstrP));                    ...
            squeeze(G(ng).x_psi(IendP,JstrP+1:JendP))';                 ...
            squeeze(flipud(G(ng).x_psi(IstrP:IendP-1,JendP)));          ...
            squeeze(fliplr(G(ng).x_psi(IstrP,JstrP:JendP-1)))'];

    Ybox = [squeeze(G(ng).y_psi(IstrP:IendP,JstrP));                    ...
            squeeze(G(ng).y_psi(IendP,JstrP+1:JendP))';                 ...
            squeeze(flipud(G(ng).y_psi(IstrP:IendP-1,JendP)));          ...
            squeeze(fliplr(G(ng).y_psi(IstrP,JstrP:JendP-1)))'];

    S.grid(ng).perimeter.X_psi = Xbox;
    S.grid(ng).perimeter.Y_psi = Ybox;

% Grid perimeter at RHO-points (counterclockwise from south). Needed for
% lateral boundary condition switch at RHO-points.

    Xbox = [squeeze(G(ng).x_rho(IstrR:IendR,JstrR));                    ...
            squeeze(G(ng).x_rho(IendR,JstrR+1:JendR))';                 ...
            squeeze(flipud(G(ng).x_rho(IstrR:IendR-1,JendR)));          ...
            squeeze(fliplr(G(ng).x_rho(IstrR,JstrR:JendR-1)))'];
    
    Ybox = [squeeze(G(ng).y_rho(IstrR:IendR,JstrR));                    ...
            squeeze(G(ng).y_rho(IendR,JstrR+1:JendR))';                 ...
            squeeze(flipud(G(ng).y_rho(IstrR:IendR-1,JendR)));          ...
            squeeze(fliplr(G(ng).y_rho(IstrR,JstrR:JendR-1)))'];

    S.grid(ng).perimeter.X_rho = Xbox;
    S.grid(ng).perimeter.Y_rho = Ybox;

% Grid perimeter at U-points (counterclockwise from south). Needed for
% lateral boundary condition switch at U-points.

    Xbox = [squeeze(G(ng).x_u(IstrU:IendU,JstrU));                      ...
            squeeze(G(ng).x_u(IendU,JstrU+1:JendU))';                   ...
            squeeze(flipud(G(ng).x_u(IstrU:IendU-1,JendU)));            ...
            squeeze(fliplr(G(ng).x_u(IstrU,JstrU:JendU-1)))'];

    Ybox = [squeeze(G(ng).y_u(IstrU:IendU,JstrU));                      ...
            squeeze(G(ng).y_u(IendU,JstrU+1:JendU))';                   ...
            squeeze(flipud(G(ng).y_u(IstrU:IendU-1,JendU)));            ...
            squeeze(fliplr(G(ng).y_u(IstrU,JstrU:JendU-1)))'];
    
    S.grid(ng).perimeter.X_u = Xbox;
    S.grid(ng).perimeter.Y_u = Ybox;

% Grid perimeter at V-points (counterclockwise from south). Needed for
% lateral boundary condition switch at V-points.

    Xbox = [squeeze(G(ng).x_v(IstrV:IendV,JstrV));                      ...
            squeeze(G(ng).x_v(IendV,JstrV+1:JendV))';                   ...
            squeeze(flipud(G(ng).x_v(IstrV:IendV-1,JendV)));            ...
            squeeze(fliplr(G(ng).x_v(IstrV,JstrV:JendV-1)))'];
    
    Ybox = [squeeze(G(ng).y_v(IstrV:IendV,JstrV));                      ...
            squeeze(G(ng).y_v(IendV,JstrV+1:JendV))';                   ...
            squeeze(flipud(G(ng).y_v(IstrV:IendV-1,JendV)));            ...
            squeeze(fliplr(G(ng).y_v(IstrV,JstrV:JendV-1)))'];

    S.grid(ng).perimeter.X_v = Xbox;
    S.grid(ng).perimeter.Y_v = Ybox;

% Grid perimeter at UV-points (counterclockwise from south). Needed for
% lateral boundary condition switch at U- and V-points. It is more
% accurate to use this perimeter with inpolygon.

    Xbox = [squeeze(G(ng).x_v(IstrV+1:IendV-1,JstrV));                  ...
            squeeze(G(ng).x_u(IendU,JstrU+1:JendU-1))';                 ...
            squeeze(flipud(G(ng).x_v(IstrV+1:IendV-1,JendV)));          ...
            squeeze(fliplr(G(ng).x_u(IstrU,JstrU+1:JendU-1)))'];
    
    Ybox = [squeeze(G(ng).y_v(IstrV+1:IendV-1,JstrV));                  ...
            squeeze(G(ng).y_u(IendU,JstrU+1:JendU-1))';                 ...
            squeeze(flipud(G(ng).y_v(IstrV+1:IendV-1,JendV)));          ...
            squeeze(fliplr(G(ng).y_u(IstrU,JstrU+1:JendU-1)))'];

    S.grid(ng).perimeter.X_uv = Xbox;
    S.grid(ng).perimeter.Y_uv = Ybox;

% Grid domain corners (PSI-points).

    Icorners = false([Im Jm]);
    Icorners(1 ,1 ) = true;
    Icorners(Im,1 ) = true;
    Icorners(Im,Jm) = true;
    Icorners(1 ,Jm) = true;
    Icorners = find(Icorners(:) == true);

    S.grid(ng).corners.index = Icorners;
    S.grid(ng).corners.X     = G(ng).x_psi(Icorners);
    S.grid(ng).corners.Y     = G(ng).y_psi(Icorners);
    S.grid(ng).corners.I     = S.grid(ng).I_psi(Icorners);
    S.grid(ng).corners.J     = S.grid(ng).J_psi(Icorners);

% Boundary edges at PSI-points, excluding corners.

    Iwest  = false([Im Jm]);    Iwest (1 ,2:Jm-1) = true;
    Isouth = false([Im Jm]);    Isouth(2:Im-1, 1) = true;
    Ieast  = false([Im Jm]);    Ieast (Im,2:Jm-1) = true;
    Inorth = false([Im Jm]);    Inorth(2:Im-1,Jm) = true;

    B(S.western_edge ).index = find(Iwest (:) == true);
    B(S.southern_edge).index = find(Isouth(:) == true);
    B(S.eastern_edge ).index = find(Ieast (:) == true);
    B(S.northern_edge).index = find(Inorth(:) == true);

    for ib=1:4,
      S.grid(ng).boundary(ib).index = B(ib).index;
      S.grid(ng).boundary(ib).X     = G(ng).x_psi(B(ib).index);
      S.grid(ng).boundary(ib).Y     = G(ng).y_psi(B(ib).index);
    end

% Boundary edges at U-points and V-points: Momentum physical boundary. 

    Iwest  = false([Im Jm+1]);  Iwest (1 ,2:Jm) = true;
    Ieast  = false([Im Jm+1]);  Ieast (Im,2:Jm) = true;

    Isouth = false([Im+1 Jm]);  Isouth(2:Im, 1) = true;
    Inorth = false([Im+1 Jm]);  Inorth(2:Im,Jm) = true;

    B(1).index = find(Iwest (:) == true);
    B(2).index = find(Isouth(:) == true);
    B(3).index = find(Ieast (:) == true);
    B(4).index = find(Inorth(:) == true);

    S.grid(ng).boundary(1).Xuv = G(ng).x_u(B(1).index);
    S.grid(ng).boundary(2).Xuv = G(ng).x_v(B(2).index);
    S.grid(ng).boundary(3).Xuv = G(ng).x_u(B(3).index);
    S.grid(ng).boundary(4).Xuv = G(ng).x_v(B(4).index);

    S.grid(ng).boundary(1).Yuv = G(ng).y_u(B(1).index);
    S.grid(ng).boundary(2).Yuv = G(ng).y_v(B(2).index);
    S.grid(ng).boundary(3).Yuv = G(ng).y_u(B(3).index);
    S.grid(ng).boundary(4).Yuv = G(ng).y_v(B(4).index);
  
  end
end


return
