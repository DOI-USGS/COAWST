function [S,G] = contact(Gnames, Cname, varargin)

%
% CONTACT:  Sets Contact Points between ROMS nested Grids.
%
% [S,G] = contact(Gnames, Cname, Lmask, MaskInterp, Lplot)
%
% This function sets contact points in the overlaping contact
% regions between nested grids. The order of nested grid file
% names in input cell array (Gnames) is important.  Set the
% file names in the order of nesting layers and time-stepping
% in ROMS.
%
% On Input:
%
%    Gnames      Input Grid NetCDF file names (cell array)
%
%    Cname       Ouptut Contact Points NetCDF file name (string)
%
%    Lmask       Switch to remove Contact Points over land
%                  (default false)
%
%    MaskInterp  Switch to interpolate PSI-, U- and V-masks (true) or
%                  computed from interpolated RHO-mask (false) using
%                  the "uvp_mask" script. We highly recommend for this
%                  switch to always be false (default false)
%
%    Lplot       Switch to plot various Contact Points figures
%                  (default false)
%
% On Output:
%
%    S           Nested grids Contact Points structure (struct array)
%    G           Nested grids structure (1 x Ngrids struct array)
%
%
% Calls to External Functions:
%
%    get_roms_grid     Gets Information Grids Structure, G(ng) 
%    grid_perimeter    Sets Nested Grids Perimeters and Boundary Edges
%    grid_connections  Sets Nested Grids Connectivity
%    plot_contact      Plots various Contact Points figures
%    write_contact     Creates and writes Contact Point data to NetCDF file
%
%
% The Contact Points structure has the following fields:
%
%    S.Ngrids                              - Number of nested grids
%    S.Ncontact                            - Number of contact regions
%    S.NLweights = 4                       - Number of linear weights
%    S.NQweights = 9                       - Number of quadratic weights
%    S.Ndatum                              - Total number of contact points
%
%    S.western_edge  = 1                   - Western  boundary edge index
%    S.southern_edge = 2                   - Southern boundary edge index
%    S.eastern_edge  = 3                   - Eastern  boundary edge index
%    S.northern_edge = 4                   - Northern boundary edge index
%
%    S.spherical                           - Spherical switch
%
%    S.grid(ng).filename                   - Grid NetCDF file name
%
%    S.grid(ng).Lp                         - Number of I-points (RHO)
%    S.grid(ng).Mp                         - Number of J-points (RHO)
%    S.grid(ng).L                          - Number of I-points (PSI)
%    S.grid(ng).M                          - Number of J-points (PSI)
%
%    S.grid(ng).refine_factor              - Refinement factor (0,3,5,7)
%
%    S.grid(ng).I_psi(:,:)                 - ROMS I-indices at PSI-points
%    S.grid(ng).J_psi(:,:)                 - ROMS J-indices at PSI-points
%    S.grid(ng).I_rho(:,:)                 - ROMS I-indices at RHO-points
%    S.grid(ng).J_rho(:,:)                 - ROMS J-indices at RHO-points
%    S.grid(ng).I_u  (:,:)                 - ROMS I-indices at U-points
%    S.grid(ng).J_u  (:,:)                 - ROMS J-indices at U-points
%    S.grid(ng).I_v  (:,:)                 - ROMS I-indices at V-points
%    S.grid(ng).J_v  (:,:)                 - ROMS I-indices at V-points
%
%    S.grid(ng).perimeter.X_psi(:)         - Perimeter X-coordinates (PSI)
%    S.grid(ng).perimeter.Y_psi(:)         - Perimeter Y-coordinates (PSI)
%    S.grid(ng).perimeter.X_rho(:)         - Perimeter X-coordinates (RHO)
%    S.grid(ng).perimeter.Y_rho(:)         - Perimeter Y-coordinates (RHO)
%    S.grid(ng).perimeter.X_u(:)           - Perimeter X-coordinates (U)
%    S.grid(ng).perimeter.Y_u(:)           - Perimeter Y-coordinates (U)
%    S.grid(ng).perimeter.X_v(:)           - Perimeter X-coordinates (V)
%    S.grid(ng).perimeter.Y_v(:)           - Perimeter Y-coordinates (V)
%                                            (counterclockwise)
%
%    S.grid(ng).corners.index(:)           - Corners linear IJ-index
%    S.grid(ng).corners.X(:)               - Corners X-coordinates (PSI)
%    S.grid(ng).corners.Y(:)               - Corners Y-coordinates (PSI)
%    S.grid(ng).corners.I(:)               - Corners I-indices (PSI)
%    S.grid(ng).corners.J(:)               - Corners J-indices (PSI)
%
%    S.grid(ng).boundary(ib).index(:)      - Boundary linear IJ-index
%    S.grid(ng).boundary(ib).X(:)          - Boundary X-coordinates (PSI)
%    S.grid(ng).boundary(ib).Y(:)          - Boundary Y-coordinates (PSI)
%                                            (without corner points)
%
%    S.contact(cr).donor_grid              - Donor grid number
%    S.contact(cr).receiver_grid           - Receiver grid number
%    S.contact(cr).coincident              - Coincident boundary switch 
%    S.contact(cr).composite               - Composite grid switch
%    S.contact(cr).mosaic                  - Mosaic grid switch
%    S.contact(cr).refinement              - Refinement grid switch
%
%    S.contact(cr).interior.okay           - true/false logical
%
%    S.contact(cr).interior.Xdg(:)         - (X,Y) coordinates and (I,J)
%    S.contact(cr).interior.Ydg(:)           indices of donor grid points
%    S.contact(cr).interior.Idg(:)           inside the receiver grid
%    S.contact(cr).interior.Jdg(:)           perimeter, [] if false 
%
%    S.contact(cr).corners.okay            - true/false logical
%
%    S.contact(cr).corners.Xdg(:)          - (X,Y) coordinates and (I,J)
%    S.contact(cr).corners.Ydg(:)            indices of donor grid points
%    S.contact(cr).corners.Idg(:)            corners laying on receiver
%    S.contact(cr).corners.Idg(:)            grid perimeter, [] if false
%
%    S.contact(cr).boundary(ib).okay       - true/false logical
%    S.contact(cr).boundary(ib).match(:)   - Donor matching points logical
%
%    S.contact(cr).boundary(ib).Xdg(:)     - (X,Y) coordinates and (I,J)
%    S.contact(cr).boundary(ib).Ydg(:)       indices of donor boundary
%    S.contact(cr).boundary(ib).Idg(:)       points laying on receiver
%    S.contact(cr).boundary(ib).Jdg(:)       grid perimeter, [] if false
%
%    S.contact(cr).point.xrg_rho(:)        - Receiver grid contact points
%    S.contact(cr).point.erg_rho(:)          (XI,ETA) coordinates, (X,Y)
%    S.contact(cr).point.Xrg_rho(:)          physical coordinates, and
%    S.contact(cr).point.Yrg_rho(:)          (I,J) indices (RHO-points)
%    S.contact(cr).point.Irg_rho(:)          where data is needed from
%    S.contact(cr).point.Jrg_rho(:)          donor
%
%    S.contact(cr).point.Idg_rho(:)        - Donor (I,J) cell containing
%    S.contact(cr).point.Jdg_rho(:)          receiver grid contact point
%
%    S.contact(cr).point.xrg_u(:)          - Receiver grid contact points
%    S.contact(cr).point.erg_u(:)            (XI,ETA) coordinates, (X,Y)
%    S.contact(cr).point.Xrg_u(:)            physical coordinates, and
%    S.contact(cr).point.Yrg_u(:)            (I,J) indices (U-points)
%    S.contact(cr).point.Irg_u(:)            where data is needed from
%    S.contact(cr).point.Jrg_u(:)            donor
%
%    S.contact(cr).point.Idg_u(:)          - Donor (I,J) cell containing
%    S.contact(cr).point.Jdg_u(:)            receiver grid contact point
%
%    S.contact(cr).point.xrg_v(:)          - Receiver grid contact points
%    S.contact(cr).point.erg_v(:)            (XI,ETA) coordinates, (X,Y)
%    S.contact(cr).point.Xrg_v(:)            physical coordinates, and
%    S.contact(cr).point.Yrg_v(:)            (I,J) indices (V-points)
%    S.contact(cr).point.Irg_v(:)            where data is needed from
%    S.contact(cr).point.Jrg_v(:)            donor
%
%    S.contact(cr).point.Idg_v(:)          - Donor (I,J) cell containing
%    S.contact(cr).point.Jdg_v(:)            receiver grid contact point
%
%    S.contact(cr).point.boundary_rho      - Contact point on RHO-boundary
%    S.contact(cr).point.boundary_u        - Contact point on   U-boundary
%    S.contact(cr).point.boundary_v        - Contact point on   V-boundary
%
%                                          - Donor data at contact point:
%    S.contact(cr).point.angle(:)              angle between XI-axis & East
%    S.contact(cr).point.f(:)                  Coriolis parameter
%    S.contact(cr).point.h(:)                  bathymetry
%    S.contact(cr).point.pm(:)                 curvilinear  XI-metric, 1/dx
%    S.contact(cr).point.pn(:)                 curvilinear ETA-metric, 1/dy
%    S.contact(cr).point.dndx(:)               d(pn)/d(xi)  metric
%    S.contact(cr).point.dmde(:)               d(pm)/d(eta) metric
%    S.contact(cr).point.mask_rho(:)           Land/Sea mask at RHO-points
%    S.contact(cr).point.mask_u(:)             Land/Sea mask at U-points
%    S.contact(cr).point.mask_v(:)             land/Sea mask at V-contact
%
%    S.refined(cr).xi_rho(:,:)             - Receiver grid curvilinear
%    S.refined(cr).eta_rho(:,:)              (XI,ETA) coordinates
%
%    S.refined(cr).x_rho(:,:)              - Receiver grid Cartesian (X,Y)
%    S.refined(cr).y_rho(:,:)                coordinates of the contact
%    S.refined(cr).x_psi(:,:)                points that requires data from
%    S.refined(cr).y_psi(:,:)                donor grid for each C-grid
%    S.refined(cr).x_u(:,:)                  type variable (RHO-, PSI, U-,
%    S.refined(cr).y_u(:,:)                  V-points). It is used to set
%    S.refined(cr).x_v(:,:)                  contact points and weights
%    S.refined(cr).y_v(:,:)                  in refinement grid.
%                                    
%    S.refined(cr).lon_rho(:,:)            - Receiver grid spherical
%    S.refined(cr).lat_rho(:,:)              (lon,lat) coordinates of the
%    S.refined(cr).lon_psi(:,:)              contact points that require
%    S.refined(cr).lat_psi(:,:)              data from donor grid for each
%    S.refined(cr).lon_u(:,:)                C-grid type variable (RHO-,
%    S.refined(cr).lat_u(:,:)                PSI, U-, and V-points). It is
%    S.refined(cr).lon_v(:,:)                to set contact points and
%    S.refined(cr).lat_v(:,:)                weights in refinement grids.
%
%    S.refined(cr).Irg_rho(:,:)            - Receiver grid ROMS indices
%    S.refined(cr).Jrg_rho(:,:)              (I,J) of the contact points
%    S.refined(cr).Irg_psi(:,:)              that require data from donor
%    S.refined(cr).Jrg_psi(:,:)              grid for each C-grid type
%    S.refined(cr).Irg_u(:,:)                variable (RHO-, PSI, U-, and 
%    S.refined(cr).Jrg_u(:,:)                V-points). It is used to set
%    S.refined(cr).Irg_v(:,:)                contact points and weights
%    S.refined(cr).Jrg_v(:,:)                in refinement grids.
%
%    S.refined(cr).mask_rho(:,:)           - Receiver grid land/sea mask
%    S.refined(cr).mask_psi(:,:)             at contact points interpolated
%    S.refined(cr).mask_u(:,:)               from donor grid at RHO-, PSI-,
%    S.refined(cr).mask_v(:,:)               U- and V-points.
%
%    S.Lweights(cr).H_rho(4,:)             - Linear weights (H) to 
%    S.Lweights(cr).H_u(4,:)                 horizontally interpolate
%    S.Lweights(cr).H_v(4,:)                 reciever data from donor grid
%
%    S.Qweights(cr).H_rho(9,:)             - Quadratic weights (H) to
%    S.Qweights(cr).H_u(9,:)                 horizontally interpolate
%    S.Qweights(cr).H_v(9,:)                 reciever data from donor grid
%
% The "refined" sub-structure is only relevant when processing the contact
% region of a refinement grid. Otherwise, it will be empty.  The setting
% of Land/Sea masking in the contact region is critical. Usually, the
% RHO-mask is interpolated from the coarser grid and the U- and V-masks
% are computed from the interpolated RHO-mask using "uvp_masks".  The 
% "MaskInterp" switch can be used to either interpolate (true) their
% values or compute using "uvp_masks" (false).  Recall that we are not
% modifying the original refined grid mask, just computing the mask in
% the contact region adjacent to the finer grid from the coarser grid mask.
% This is only relevant when there are land/sea masking features in any
% of the refinement grid physical boundaries. If so, the user just needs to
% experiment with "MaskInterp" and edit such points during post-processing.
%
% The locations of the spatial linear interpolation weights in the donor
% grid with respect the receiver grid contact region at a particulat
% contact point x(Irg,Jrg,Krg) are:
%
%                       8___________7   (Idg+1,Jdg+1,Kdg)
%                      /.          /|
%                     / .         / |
%  (Idg,Jdg+1,Kdg)  5/___________/6 |
%                    |  .        |  |
%                    |  .   x    |  |
%                    | 4.........|..|3  (Idg+1,Jdg+1,Kdg-1)
%                    | .         |  /
%                    |.          | /
%                    |___________|/
%  (Idg,Jdg,Kdg-1)   1           2
%
%                                        Suffix:   dg = donor grid
%                                                  rg = receiver grid
%
% We just need to set the horizontal interpolation weights 1:4. The
% Other weights needed for vertical interpolation will be set inside
% ROMS.  Notice that if the contact point "cp" between donor and
% receiver grids are coincident:
%
%          S.Lweights(cr).H_rho(1,cp) = 1.0 
%          S.Lweights(cr).H_rho(2,cp) = 0.0 
%          S.Lweights(cr).H_rho(3,cp) = 0.0 
%          S.Lweights(cr).H_rho(4,cp) = 0.0 
% Then
%          receiver_value(Irg,Jrg) = donor_value(Idg,Jdg)
%

% svn $Id: contact.m 746 2014-12-15 23:24:27Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

Lmask = false;
MaskInterp = false;
Lplot = false;

switch numel(varargin)
  case 1
    Lmask = varargin{1};
  case 2
    Lmask = varargin{1};
    MaskInterp = varargin{2};
  case 3
    Lmask = varargin{1};
    MaskInterp = varargin{2};
    Lplot = varargin{3};
end
  
Ngrids = length(Gnames);           % number of nested grids
Ncontact = (Ngrids-1)*2;           % number of contact regions

%--------------------------------------------------------------------------
% Get nested grids information and set perimeters and boundary edges.
%--------------------------------------------------------------------------

% Get nested grid structures.

G = grids_structure(Gnames);

% Set nested grids perimeters and boundary edges.

S = grid_perimeter(G);

%--------------------------------------------------------------------------
% Set up nested grid connections.
%--------------------------------------------------------------------------

S = grid_connections(G, S);

%--------------------------------------------------------------------------
% Set contact points in each contact region.
%--------------------------------------------------------------------------

Ndatum = 0;

for cr=1:Ncontact,
  dg = S.contact(cr).donor_grid;
  rg = S.contact(cr).receiver_grid;
  if (S.contact(cr).coincident),
    P = coincident (cr, dg, rg, Lmask, G, S);
    S.contact(cr).point = P;
  elseif (S.contact(cr).refinement),
    [P, R] = refinement (cr, dg, rg, Lmask, G, S, MaskInterp);
    S.contact(cr).point = P;
    S.refined(cr) = R;
  end

  Ndatum = Ndatum + length(S.contact(cr).point.Irg_rho)                 ...
                  + length(S.contact(cr).point.Irg_u)                   ...
                  + length(S.contact(cr).point.Irg_v);
end

S.Ndatum = Ndatum;

%--------------------------------------------------------------------------
% Determine which contact points lay on the receiver grid boundary. This
% information is used in ROMS to avoid applying lateral boundary
% conditions at those points.
%--------------------------------------------------------------------------

S = boundary_contact(S);

%--------------------------------------------------------------------------
% Set contact points horizontal interpolation weights.
%--------------------------------------------------------------------------
%
% The weights are not scaled according to the land/sea mask.  This is done
% in ROMS because the land/sea masking time dependecy when wetting and
% drying is activated for an application.

ImposeMask = false;

S = Hweights(G, S, ImposeMask);

%--------------------------------------------------------------------------
% Create and write out Contact Point data into output NetCDF file.
%--------------------------------------------------------------------------

write_contact(Cname, S, G);

%--------------------------------------------------------------------------
% Plot contact areas and contact points.
%--------------------------------------------------------------------------

if (Lplot),
  plot_contact(G, S);
end

return

%--------------------------------------------------------------------------

function R = refine_coordinates(cr, dg, rg, G, S, MaskInterp)

%
% R = refine_coordinates(cr, dg, rg, G, S, MaskInterp)
%
% This function computes receiver grid refinement coordinates from the
% donor grid.  The refinement coordinates stencil is larger to facilitate
% the extraction of receiver grid contact points coordinates elsewhere.
%
% On Input:
%
%    cr          Contact region number
%    dg          Donor grid number
%    rg          Receiver grid number
%    G           Nested grids structure (1 x Ngrid struct array)
%    S           Nested grids information structure (struct array)
%    MaskInterp  Switch to interpolate PSI-, U- and V-mask (true) or
%                  computed from interpolated RHO-mask (false)
% On Output:
%
%    R           Refinement coordinates structure (struct array)
%
%                  R.xi_rho         - Receiver grid curvilinear (XI,ETA)
%                  R.eta_rho          coordinates in terms of the donor
%                  R.xi_psi           grid coordinates for each C-grid
%                  R.eta_psi          type variable (RHO-, PSI-, U-, and
%                  R.xi_u             V-points). It will be used elsewhere
%                  R.eta_u            for setting contact points and
%                  R.xi_v             interpolation weights more
%                  R.eta_v            accuratelly.
%
%                if (~spherical)
%
%                  R.x_rho          - Receiver grid Cartesian (X,Y)
%                  R.y_rho            coordinates of the contact
%                  R x_psi            points that requires data from
%                  R y_psi            donor grid for each C-grid type
%                  R.x_u              variable (RHO-, PSI-, U-, and 
%                  R.y_u              V-points). It will be used
%                  R.x_v              elsewhere for setting contact
%                  R.y_v              points and interpolation weights.
%
%                elseif (spherical)
%                                    
%                  R.lon_rho        - Receiver grid spherical (lon,lat)
%                  R.lat_rho          coordinates of the contact
%                  R lon_psi          points that requires data from
%                  R lat_psi          donor grid for each C-grid type
%                  R.lon_u            variable (RHO-, PSI-, U-, and 
%                  R.lat_u            V-points). It will be used
%                  R.lon_v            elsewhere for setting contact
%                  R.lat_v            points and interpolation weights.
%
%                end
%
%                  R.Irg_rho        - Receiver grid (I,J) coordinates
%                  R.Jrg_rho          (ROMS indices) of the contact
%                  R Irg_psi          points that requires data from
%                  R Jrg_psi          donor grid for each C-grid type
%                  R.Irg_u            variable (RHO-, PSI-, U-, and 
%                  R.Jrg_u            V-points). It will be used
%                  R.Irg_v            elsewhere for setting contact
%                  R.Jrg_v            points and interpolation weights.
%
%                  R.mask_rho       - Receiver grid land/sea mask of
%                  R.mask_psi         the contact points interpolated
%                  R.mask_u           from donor grid at RHO-, PSI-,
%                  R.mask_v           U- and V-points.
%  
% Notice the land/sea mask of the contact points interpolated from the
% donor grid may be require some manual editing during post-processing
% of contact points data.
%

% Initialize.

spherical = S.spherical;

% Set TriScatteredInterp method.  Use 'natural' for Natural Neighbor
% Interpolation Method.  This does NOT implies Nearest Neighbor
% Interpolation.  Natural Neighbor Interpolation is a triangulantion
% based method that has an area-of-influence weighting associated
% with each sample data point.  It is C1 continous except at the
% sample locations.

method = 'natural';

%--------------------------------------------------------------------------
% Set extraction indices with respect the donor (coarse) grid.
%--------------------------------------------------------------------------

% Get donor and reciever grid dimension parameters.

Lp = S.grid(dg).Lp;   Mp = S.grid(dg).Mp;
L  = S.grid(dg).L;    M  = S.grid(dg).M;

% Set donor (coarse) grid fractional coordinates.

[YrC, XrC] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);      % RHO-points
[YpC, XpC] = meshgrid(1.0:1:M     , 1.0:1:L     );      % PSI-points
[YuC, XuC] = meshgrid(0.5:1:Mp-0.5, 1.0:1:L     );      % U-points
[YvC, XvC] = meshgrid(1.0:1:M     , 0.5:1:Lp-0.5);      % V-points

% Set extraction fractional coordinates from donor grid.  The corners
% values in the connectivity structure, C, for the current contact
% region are actually the corners of the refinement  grid physical
% perimeter (PSI-points). Therefore, we need to extract a larger grid
% in terms of the following offsets:
%
%   3 extra PSI-points adjacent to the western/southern perimeter edges
%   2 extra PSI-points adjacent to the eastern/northern perimeter edges
%
% The other C-grid types variables will be accounted by the "fac0" and
% "fac1" factors.

delta = 1.0/S.grid(rg).refine_factor;
half  = 0.5*delta;

offset_west  = 3.0*delta;
offset_south = 3.0*delta;
offset_east  = 2.0*delta;
offset_north = 2.0*delta;

Imin = min(S.contact(cr).corners.Idg)-offset_west;
Imax = max(S.contact(cr).corners.Idg)+offset_east;
Jmin = min(S.contact(cr).corners.Jdg)-offset_south;
Jmax = max(S.contact(cr).corners.Jdg)+offset_north;
  
% Set receiver (fine) grid fractional coordinates.

IpF = (Imin:delta:Imax);                                    % PSI-points
JpF = (Jmin:delta:Jmax);                                    % PSI-points
IrF = [IpF(1)-half IpF+half];                               % RHO-points
JrF = [JpF(1)-half JpF+half];                               % RHO-points

[YrF, XrF] = meshgrid(JrF, IrF);                            % RHO-points
[YpF, XpF] = meshgrid(JpF, IpF);                            % PSI-points
[YuF, XuF] = meshgrid(JrF, IpF);                            % U-points
[YvF, XvF] = meshgrid(JpF, IrF);                            % V-points

%--------------------------------------------------------------------------
% Extract larger receiver (fine) grid from donor (coarse) grid.
%--------------------------------------------------------------------------

R.spherical = spherical;
R.uniform   = G(dg).uniform;

% Set (XI,ETA) coordinates for the receiver finer grid in terms of the
% donor coarser grid.  This will be needed for more precise determination
% of the interpolation weights.

% RHO-points coordinates.

R.xi_rho  = XrF;
R.eta_rho = YrF;

if (spherical),
  if (~isempty(G(dg).x_rho) && ~isempty(G(dg).y_rho)),

    FCr = TriScatteredInterp(XrC(:), YrC(:),                            ...
                             G(dg).x_rho(:), method);

                               R.x_rho = FCr(XrF, YrF);
    FCr.V = G(dg).y_rho(:);    R.y_rho = FCr(XrF, YrF); 

    FSr = TriScatteredInterp(G(dg).x_rho(:), G(dg).y_rho(:),            ...
                             G(dg).lon_rho(:), method);

                               R.lon_rho = FSr(R.x_rho, R.y_rho);
    FSr.V = G(dg).lat_rho(:);  R.lat_rho = FSr(R.x_rho, R.y_rho);
   
  elseif (~isempty(G(dg).lon_rho) && ~isempty(G(dg).lat_rho)),

    FSr = TriScatteredInterp(XrC(:), YrC(:),                            ...
                             G(dg).lon_rho(:), method);

                               R.lon_rho = FSr(XrF, YrF);
    FSr.V = G(dg).lat_rho(:);  R.lat_rho = FSr(XrF, YrF);

    FCr = TriScatteredInterp(G(dg).lon_rho(:), G(dg).lat_rho(:),        ...
                             G(dg).x_rho(:),method);

                               R.x_rho = RCr(R.lon_rho, R.lat_rho);
    FCr.V = G(dg).y_rho(:);    R.y_rho = RCr(R.lon_rho, R.lat_rho);  
  end
  
  [Im, Jm] = size(R.lon_rho);
  [Jrg, Irg] = meshgrid(1:Jm, 1:Im);
  R.Irg_rho = Irg - 4;                              % ROMS indices:
  R.Jrg_rho = Jrg - 4;                              % -3:L+2, -3:M+2
  
else

  FCr = TriScatteredInterp(XrC(:), YrC(:),                              ...
                           G(dg).x_rho(:), method);

                               R.x_rho = FCr(XrF, YrF);
  FCr.V = G(dg).y_rho(:);      R.y_rho = FCr(XrF, YrF); 

  [Im, Jm] = size(R.x_rho);
  [Jrg, Irg] = meshgrid(1:Jm, 1:Im);
  R.Irg_rho = Irg - 4;                               % ROMS indices:
  R.Jrg_rho = Jrg - 4;                               % -3:L+2, -3:M+2

end

% PSI-points coordinates.

R.xi_psi  = XpF;
R.eta_psi = YpF;

if (spherical),
  if (~isempty(G(dg).x_psi) && ~isempty(G(dg).y_psi)),

    FCp = TriScatteredInterp(XpC(:), YpC(:),                            ...
                             G(dg).x_psi(:), method);

                               R.x_psi = FCp(XpF, YpF);
    FCp.V = G(dg).y_psi(:);    R.y_psi = FCp(XpF, YpF); 

    FSp = TriScatteredInterp(G(dg).x_psi(:), G(dg).y_psi(:),            ...
                             G(dg).lon_psi(:), method);

                               R.lon_psi = FSp(R.x_psi, R.y_psi);
    FSp.V = G(dg).lat_psi(:);  R.lat_psi = FSp(R.x_psi, R.y_psi);
  
  elseif (~isempty(G(dg).lon_psi) && ~isempty(G(dg).lat_psi)),

    FSp = TriScatteredInterp(XpC(:), YpC(:),                            ...
                             G(dg).lon_psi(:), method);

                               R.lon_psi = FSp(XpF, YpF);
    FSp.V = G(dg).lat_psi(:);  R.lat_psi = FSp(XpF, YpF);

    FCp = TriScatteredInterp(G(dg).lon_psi(:), G(dg).lat_psi(:),        ...
                             G(dg).x_psi(:), method);

                               R.x_rho = RCp(R.lon_rho, R.lat_rho);
    FCp.V = G(dg).y_rho(:);    R.y_rho = RCp(R.lon_rho, R.lat_rho);  
  end

  [Im, Jm] = size(R.lon_psi);
  [Jrg, Irg] = meshgrid(1:Jm, 1:Im);
  R.Irg_psi = Irg - 3;                               % ROMS indices:
  R.Jrg_psi = Jrg - 3;                               % -2:L+2, -2:M+2

else

  FCp = TriScatteredInterp(XpC(:), YpC(:),                              ...
                           G(dg).x_psi(:), method);

                               R.x_psi = FCp(XpF, YpF);
  FCp.V = G(dg).y_psi(:);      R.y_psi = FCp(XpF, YpF); 

  [Im, Jm] = size(R.x_psi);
  [Jrg, Irg] = meshgrid(1:Jm, 1:Im);
  R.Irg_psi = Irg - 3;                               % ROMS indices:
  R.Jrg_psi = Jrg - 3;                               % -2:L+2, -2:M+2
end

% U-points coordinates.

R.xi_u  = XuF;
R.eta_u = YuF;

if (spherical),
  if (~isempty(G(dg).x_u) && ~isempty(G(dg).y_u)),

    FCu = TriScatteredInterp(XuC(:), YuC(:),                            ...
                             G(dg).x_u(:), method);

                               R.x_u = FCu(XuF, YuF);
    FCu.V = G(dg).y_u(:);      R.y_u = FCu(XuF, YuF); 

    FSu = TriScatteredInterp(G(dg).x_u(:), G(dg).y_u(:),                ...
                             G(dg).lon_u(:), method);

                               R.lon_u = FSu(R.x_u, R.y_u);
    FSu.V = G(dg).lat_u(:);    R.lat_u = FSu(R.x_u, R.y_u);
  
  elseif (~isempty(G(dg).lon_u) && ~isempty(G(dg).lat_u)),

    FSu = TriScatteredInterp(XuC(:), YuC(:),                            ...
                             G(dg).lon_u(:), method);

                               R.lon_u = FSu(XuF, YuF);
    FSu.V = G(dg).lat_u(:);    R.lat_u = FSu(XuF, YuF);

    FCu = TriScatteredInterp(G(dg).lon_u(:), G(dg).lat_u(:),            ...
                             G(dg).x_u(:), method);

                               R.x_u = FCu(R.lon_u, R.lat_u);
    FCu.V = G(dg).y_u(:);      R.y_u = FCu(R.lon_u, R.lat_u);  
  end

  [Im, Jm]=size(R.lon_u);
  [Jrg, Irg] = meshgrid(1:Jm, 1:Im);
  R.Irg_u = Irg - 3;                                % ROMS indices:
  R.Jrg_u = Jrg - 4;                                % -2:L+2, -3:M+2

else

  FCu = TriScatteredInterp(XuC(:), YuC(:),                              ...
                           G(dg).x_u(:), method);

                              R.x_u = FCu(XuF, YuF);
  FCu.V = G(dg).y_u(:);       R.y_u = FCu(XuF, YuF); 

  [Im, Jm] = size(R.x_u);
  [Jrg, Irg] = meshgrid(1:Jm, 1:Im);
  R.Irg_u = Irg - 3;                                % ROMS indices:
  R.Jrg_u = Jrg - 4;                                % -2:L+2, -3:M+2

end

% V-points coordinates.

R.xi_v  = XvF;
R.eta_v = YvF;

if (spherical),
  if (~isempty(G(dg).x_v) && ~isempty(G(dg).y_v)),

    FCv = TriScatteredInterp(XvC(:), YvC(:),                            ...
                             G(dg).x_v(:), method);

                               R.x_v = FCv(XvF, YvF);
    FCv.V = G(dg).y_v(:);      R.y_v = FCv(XvF, YvF); 

    FSv = TriScatteredInterp(G(dg).x_v(:), G(dg).y_v(:),                ...
                             G(dg).lon_v(:), method);

                               R.lon_v = FSv(R.x_v, R.y_v);
    FSv.V = G(dg).lat_v(:);    R.lat_v = FSv(R.x_v, R.y_v);
  
  elseif (~isempty(G(dg).lon_v) && ~isempty(G(dg).lat_v)),

    FSv = TriScatteredInterp(XvC(:), YvC(:),                            ...
                             G(dg).lon_v(:), method);

                               R.lon_v = FSv(XvF, YvF);
    FSv.V = G(dg).lat_v(:);    R.lat_v = FSv(XvF, YvF);

    FCv = TriScatteredInterp(G(dg).lon_v(:), G(dg).lat_v(:),            ...
                             G(dg).x_v(:), method);

                               R.x_v = FCv(R.lon_v, R.lat_v);
    FCv.V = G(dg).y_v(:);      R.y_v = FCv(R.lon_v, R.lat_v);  
  end

  [Im, Jm] = size(R.lon_v);
  [Jrg, Irg] = meshgrid(1:Jm, 1:Im);
  R.Irg_v = Irg - 4;                                % ROMS indices:
  R.Jrg_v = Jrg - 3;                                % -3:L+2, -2:M+2

else

  FCv = TriScatteredInterp(XvC(:), YvC(:),                              ...
                           G(dg).x_v(:), method);

                               R.x_v = FCv(XvF, YvF);
  FCv.V = G(dg).y_v(:);        R.y_v = FCv(XvF, YvF); 

  [Im, Jm] = size(R.x_v);
  [Jrg, Irg] = meshgrid(1:Jm, 1:Im);
  R.Irg_v = Irg - 4;                                % ROMS indices:
  R.Jrg_v = Jrg - 3;                                % -3:L+2, -2:M+2

end

% Interpolate other grid variables.

if (spherical),
  if (~isempty(G(dg).x_rho) && ~isempty(G(dg).y_rho)),

    FCr.V = G(dg).angle(:);    R.angle = FCr(XrF, YrF); 
    FCr.V = G(dg).f(:);        R.f     = FCr(XrF, YrF); 
    FCr.V = G(dg).h(:);        R.h     = FCr(XrF, YrF); 

  elseif (~isempty(G(dg).lon_rho) && ~isempty(G(dg).lat_rho)),

    FSr.V = G(dg).angle(:);    R.angle = FCr(XrF, YrF); 
    FSr.V = G(dg).f(:);        R.f     = FCr(XrF, YrF); 
    FSr.V = G(dg).h(:);        R.h     = FCr(XrF, YrF); 
  
  end

else

  FCr.V = G(dg).angle(:);      R.angle = FCr(XrF, YrF); 
  FCr.V = G(dg).f(:);          R.f     = FCr(XrF, YrF); 
  FCr.V = G(dg).h(:);          R.h     = FCr(XrF, YrF); 

end  

% Set grid metrics.

if (G(dg).uniform),

% Donor has a uniform grid distribution.

  dx = 1/unique(G(dg).pm(:));
  dy = 1/unique(G(dg).pn(:));

  R.pm = ones(size(R.h)) .* (S.grid(rg).refine_factor/dx);
  R.pn = ones(size(R.h)) .* (S.grid(rg).refine_factor/dy);
  
  R.dndx = zeros(size(R.h));
  R.dmde = zeros(size(R.h));

else

% Recompute metrics at refined resolution.  We cannot interpolate.
% The computation of metrics from discrete point is subject to
% roundoff.  There is not much that we can do here.  The roundoff
% is small and of the order 1.0E-16 (eps value).

  if (spherical),
    GreatCircle = true;
  else
    GreatCircle = false;
  end

  [R.pm, R.pn, R.dndx, R.dmde] = grid_metrics(R, GreatCircle);

end

% Get Land/sea masking. In some applications is better to linearly
% interpolate the finer grid mask because recomputing the mask in
% the contact region may give different values. Sometimes, the
% U- and V-mask is one point smaller or larger. We are not
% modifying the original refined grid mask, but just computing
% the mask in the contact region adjacent to finer grid from the
% coarser grid mask. The user just need to experiment with this
% and edit such points during post-processing.

R.mask_rho = interp2(XrC', YrC', G(dg).mask_rho', XrF, YrF, 'nearest');

if (MaskInterp),
  R.mask_psi = interp2(XpC', YpC', G(dg).mask_psi', XpF, YpF, 'nearest');
  R.mask_u   = interp2(XuC', YuC', G(dg).mask_u'  , XuF, YuF, 'nearest');
  R.mask_v   = interp2(XvC', YvC', G(dg).mask_v'  , XvF, YvF, 'nearest');
else
  [umask, vmask, pmask] = uvp_masks(R.mask_rho);
  R.mask_psi = pmask;
  R.mask_u   = umask;
  R.mask_v   = vmask;
end  

% Over-write boundary values with original data.  This is done to avoid
% roundoff in "inpolygon" so we can match boundary points exactly.

Lp = S.grid(rg).Lp;   Mp = S.grid(rg).Mp;
L  = S.grid(rg).L;    M  = S.grid(rg).M;

IstrR = 1;      IstrP = 1;     IstrU = 1;      IstrV = 1;
IendR = Lp;     IendP = L;     IendU = L;      IendV = Lp;
JstrR = 1;      JstrP = 1;     JstrU = 1;      JstrV = 1;
JendR = Mp;     JendP = M;     JendU = Mp;     JendV = M;

if (spherical),
  R.lon_psi(IstrP+3:IendP+3,JstrP+3)=G(rg).lon_psi(IstrP:IendP,JstrP);
  R.lon_psi(IstrP+3:IendP+3,JendP+3)=G(rg).lon_psi(IstrP:IendP,JendP);
  R.lon_psi(IstrP+3,JstrP+4:JendP+2)=G(rg).lon_psi(IstrP,JstrP+1:JendP-1);
  R.lon_psi(IendP+3,JstrP+4:JendP+2)=G(rg).lon_psi(IendP,JstrP+1:JendP-1);

  R.lat_psi(IstrP+3:IendP+3,JstrP+3)=G(rg).lat_psi(IstrP:IendP,JstrP);
  R.lat_psi(IstrP+3:IendP+3,JendP+3)=G(rg).lat_psi(IstrP:IendP,JendP);
  R.lat_psi(IstrP+3,JstrP+4:JendP+2)=G(rg).lat_psi(IstrP,JstrP+1:JendP-1);
  R.lat_psi(IendP+3,JstrP+4:JendP+2)=G(rg).lat_psi(IendP,JstrP+1:JendP-1);

  R.lon_rho(IstrR+3:IendR+3,JstrR+3)=G(rg).lon_rho(IstrR:IendR,JstrR);
  R.lon_rho(IstrR+3:IendR+3,JendR+3)=G(rg).lon_rho(IstrR:IendR,JendR);
  R.lon_rho(IstrR+3,JstrR+4:JendR+2)=G(rg).lon_rho(IstrR,JstrR+1:JendR-1);
  R.lon_rho(IendR+3,JstrR+4:JendR+2)=G(rg).lon_rho(IendR,JstrR+1:JendR-1);

  R.lat_rho(IstrR+3:IendR+3,JstrR+3)=G(rg).lat_rho(IstrR:IendR,JstrR);
  R.lat_rho(IstrR+3:IendR+3,JendR+3)=G(rg).lat_rho(IstrR:IendR,JendR);
  R.lat_rho(IstrR+3,JstrR+4:JendR+2)=G(rg).lat_rho(IstrR,JstrR+1:JendR-1);
  R.lat_rho(IendR+3,JstrR+4:JendR+2)=G(rg).lat_rho(IendR,JstrR+1:JendR-1);

  R.lon_u(IstrU+3:IendU+3,JstrU+3)=G(rg).lon_u(IstrU:IendU,JstrU);
  R.lon_u(IstrU+3:IendU+3,JendU+3)=G(rg).lon_u(IstrU:IendU,JendU);
  R.lon_u(IstrU+3,JstrU+4:JendU+2)=G(rg).lon_u(IstrU,JstrU+1:JendU-1);
  R.lon_u(IendU+3,JstrU+4:JendU+2)=G(rg).lon_u(IendU,JstrU+1:JendU-1);
 
  R.lat_u(IstrU+3:IendU+3,JstrU+3)=G(rg).lat_u(IstrU:IendU,JstrU);
  R.lat_u(IstrU+3:IendU+3,JendU+3)=G(rg).lat_u(IstrU:IendU,JendU);
  R.lat_u(IstrU+3,JstrU+4:JendU+2)=G(rg).lat_u(IstrU,JstrU+1:JendU-1);
  R.lat_u(IendU+3,JstrU+4:JendU+2)=G(rg).lat_u(IendU,JstrU+1:JendU-1);

  R.lon_v(IstrV+3:IendV+3,JstrV+3)=G(rg).lon_v(IstrV:IendV,JstrV);
  R.lon_v(IstrV+3:IendV+3,JendV+3)=G(rg).lon_v(IstrV:IendV,JendV);
  R.lon_v(IstrV+3,JstrV+4:JendV+2)=G(rg).lon_v(IstrV,JstrV+1:JendV-1);
  R.lon_v(IendV+3,JstrV+4:JendV+2)=G(rg).lon_v(IendV,JstrV+1:JendV-1);

  R.lat_v(IstrV+3:IendV+3,JstrV+3)=G(rg).lat_v(IstrV:IendV,JstrV);
  R.lat_v(IstrV+3:IendV+3,JendV+3)=G(rg).lat_v(IstrV:IendV,JendV);
  R.lat_v(IstrV+3,JstrV+4:JendV+2)=G(rg).lat_v(IstrV,JstrV+1:JendV-1);
  R.lat_v(IendV+3,JstrV+4:JendV+2)=G(rg).lat_v(IendV,JstrV+1:JendV-1);

else

  R.x_psi(IstrP+3:IendP+3,JstrP+3)=G(rg).x_psi(IstrP:IendP,JstrP);
  R.x_psi(IstrP+3:IendP+3,JendP+3)=G(rg).x_psi(IstrP:IendP,JendP);
  R.x_psi(IstrP+3,JstrP+4:JendP+2)=G(rg).x_psi(IstrP,JstrP+1:JendP-1);
  R.x_psi(IendP+3,JstrP+4:JendP+2)=G(rg).x_psi(IendP,JstrP+1:JendP-1);

  R.y_psi(IstrP+3:IendP+3,JstrP+3)=G(rg).y_psi(IstrP:IendP,JstrP);
  R.y_psi(IstrP+3:IendP+3,JendP+3)=G(rg).y_psi(IstrP:IendP,JendP);
  R.y_psi(IstrP+3,JstrP+4:JendP+2)=G(rg).y_psi(IstrP,JstrP+1:JendP-1);
  R.y_psi(IendP+3,JstrP+4:JendP+2)=G(rg).y_psi(IendP,JstrP+1:JendP-1);

  R.x_rho(IstrR+3:IendR+3,JstrR+3)=G(rg).x_rho(IstrR:IendR,JstrR);
  R.x_rho(IstrR+3:IendR+3,JendR+3)=G(rg).x_rho(IstrR:IendR,JendR);
  R.x_rho(IstrR+3,JstrR+4:JendR+2)=G(rg).x_rho(IstrR,JstrR+1:JendR-1);
  R.x_rho(IendR+3,JstrR+4:JendR+2)=G(rg).x_rho(IendR,JstrR+1:JendR-1);

  R.y_rho(IstrR+3:IendR+3,JstrR+3)=G(rg).y_rho(IstrR:IendR,JstrR);
  R.y_rho(IstrR+3:IendR+3,JendR+3)=G(rg).y_rho(IstrR:IendR,JendR);
  R.y_rho(IstrR+3,JstrR+4:JendR+2)=G(rg).y_rho(IstrR,JstrR+1:JendR-1);
  R.y_rho(IendR+3,JstrR+4:JendR+2)=G(rg).y_rho(IendR,JstrR+1:JendR-1);
  
  R.x_u(IstrU+3:IendU+3,JstrU+3)=G(rg).x_u(IstrU:IendU,JstrU);
  R.x_u(IstrU+3:IendU+3,JendU+3)=G(rg).x_u(IstrU:IendU,JendU);
  R.x_u(IstrU+3,JstrU+4:JendU+2)=G(rg).x_u(IstrU,JstrU+1:JendU-1);
  R.x_u(IendU+3,JstrU+4:JendU+2)=G(rg).x_u(IendU,JstrU+1:JendU-1);
 
  R.y_u(IstrU+3:IendU+3,JstrU+3)=G(rg).y_u(IstrU:IendU,JstrU);
  R.y_u(IstrU+3:IendU+3,JendU+3)=G(rg).y_u(IstrU:IendU,JendU);
  R.y_u(IstrU+3,JstrU+4:JendU+2)=G(rg).y_u(IstrU,JstrU+1:JendU-1);
  R.y_u(IendU+3,JstrU+4:JendU+2)=G(rg).y_u(IendU,JstrU+1:JendU-1);

  R.x_v(IstrV+3:IendV+3,JstrV+3)=G(rg).x_v(IstrV:IendV,JstrV);
  R.x_v(IstrV+3:IendV+3,JendV+3)=G(rg).x_v(IstrV:IendV,JendV);
  R.x_v(IstrV+3,JstrV+4:JendV+2)=G(rg).x_v(IstrV,JstrV+1:JendV-1);
  R.x_v(IendV+3,JstrV+4:JendV+2)=G(rg).x_v(IendV,JstrV+1:JendV-1);

  R.y_v(IstrV+3:IendV+3,JstrV+3)=G(rg).y_v(IstrV:IendV,JstrV);
  R.y_v(IstrV+3:IendV+3,JendV+3)=G(rg).y_v(IstrV:IendV,JendV);
  R.y_v(IstrV+3,JstrV+4:JendV+2)=G(rg).y_v(IstrV,JstrV+1:JendV-1);
  R.y_v(IendV+3,JstrV+4:JendV+2)=G(rg).y_v(IendV,JstrV+1:JendV-1);
end

return

%--------------------------------------------------------------------------

function C = coincident(cr, dg, rg, Lmask, G, S)

%
% C = coincident(cr, dg, rg, Lmask, G, S)
%
% This function sets the contact points in the receiver grid contact
% region the requires donor grid data.  It computes the contact points
% for the case that receiver and donor grids contact points are fully
% coincident (mosaic class) or have some coincident points (composite
% class).
%
% On Input:
%
%    cr          Contact region number
%    dg          Donor grid number
%    rg          Receiver grid number
%    Lmask       Switch to remove contact points on land
%    G           Nested grids structure (1 x Ngrid struct array)
%    S           Nested grids information structure (struct array)
%
% On Output:
%
%    C           Contact points structure (struct array):
%
%                  C.xrg_rho(:)            - Receiver grid contact points
%                  C.erg_rho(:)              (XI,ETA) coordinates, (X,Y)
%                  C.Xrg_rho(:)              physical coordinates, and
%                  C.Yrg_rho(:)              and (I,J) indices (RHO-points)
%                  C.Irg_rho(:)              where data is needed from
%                  C.Jrg_rho(:)              donor
%
%                  C.Idg_rho(:)            - Donor (I,J) cell containing
%                  C.Jdg_rho(:)              receiver grid contact point
%
%                  C.xrg_u(:)              - Receiver grid contact points
%                  C.erg_u(:)                (XI,ETA) coordinates, (X,Y)
%                  C.Xrg_u(:)                physical coordinates, and
%                  C.Yrg_u(:)                and (I,J) indices (U-points)
%                  C.Irg_u(:)                where data is needed from
%                  C.Jrg_u(:)                donor
%
%                  C.Idg_u(:)              - Donor (I,J) cell containing
%                  C.Jdg_u(:)                receiver grid contact point
%
%                  C.xrg_v(:)              - Receiver grid contact points
%                  C.erg_v(:)                (XI,ETA) coordinates, (X,Y)
%                  C.Xrg_v(:)                physical coordinates, and
%                  C.Yrg_v(:)                and (I,J) indices (V-points)
%                  C.Irg_v(:)                where data is needed from
%                  C.Jrg_v(:)                donor
%
%                  C.Idg_v(:)              - Donor (I,J) cell containing
%                  C.Jdg_v(:)                receiver grid contact point
%
%                  C.boundary_rho          - Contact point on RHO-boundary
%                  C.boundary_u            - Contact point on   U-boundary
%                  C.boundary_v            - Contact point on   V-boundary
%                                            (initialized to empty)
%
%                                          - Donor data at contact point:
%                  C.angle(:)                  angle between XI-axis & East
%                  C.f(:)                      Coriolis parameter
%                  C.h(:)                      bathymetry
%                  C.pm(:)                     curvilinear  XI-metric, 1/dx
%                  C.pn(:)                     curvilinear ETA-metric, 1/dy
%                  C.dndx(:)                   d(pn)/d(xi)  metric
%                  C.dmde(:)                   d(pm)/d(eta) metric
%                  C.mask_rho(:)               Land/Sea mask at RHO-points
%                  C.mask_u(:)                 Land/Sea mask at U-points
%                  C.mask_v(:)                 land/Sea mask at V-contact
%

% Initialize.

iwest  = S.western_edge;
isouth = S.southern_edge;
ieast  = S.eastern_edge;
inorth = S.northern_edge;

Lpdg = S.grid(dg).Lp;   Mpdg = S.grid(dg).Mp;
Lprg = S.grid(rg).Lp;   Mprg = S.grid(rg).Mp;

Ldg  = S.grid(dg).L;    Mdg  = S.grid(dg).M;
Lrg  = S.grid(rg).L;    Mrg  = S.grid(rg).M;

spherical = S.spherical;

C = struct('xrg_rho'     , [], 'erg_rho'     , [],                      ...
           'Xrg_rho'     , [], 'Yrg_rho'     , [],                      ...
           'Irg_rho'     , [], 'Jrg_rho'     , [],                      ...
           'Idg_rho'     , [], 'Jdg_rho'     , [],                      ...
           'Xrg_u'       , [], 'Yrg_u'       , [],                      ...
           'Irg_u'       , [], 'Jrg_u'       , [],                      ...
           'Idg_u'       , [], 'Jdg_u'       , [],                      ...
           'Xrg_v'       , [], 'Yrg_v'       , [],                      ...
           'Irg_v'       , [], 'Jrg_v'       , [],                      ...
           'Idg_v'       , [], 'Jdg_v'       , [],                      ...
           'boundary_rho', [], 'boundary_u'  , [], 'boundary_v'  , [],  ...
           'angle'       , [], 'f'           , [], 'h'           , [],  ...
           'pm'          , [], 'pn'          , [],                      ...
           'dndx'        , [], 'dmde'        , [],                      ...
           'mask_rho'    , [], 'mask_u'      , [], 'mask_v'      , []);

W = C;

%--------------------------------------------------------------------------
% Set receiver grid contact points using the donor grid data for the
% current contact region.
%--------------------------------------------------------------------------

% Western boundary of the donor grid lies on the receiver grid eastern
% boundary.

if (S.contact(cr).boundary(iwest).okay),
  Ioffr = 3;
  Ioffu = 3;
  Ioffv = 3;

  ib = iwest;

  Ir = S.grid(dg).I_rho <= Ioffr;  Ir(1,:) = false;
  Iu = S.grid(dg).I_u   <= Ioffu;
  Iv = S.grid(dg).I_v   <= Ioffv;  Iv(1,:) = false;

  W(ib).xrg_rho = S.grid(dg).XI_rho(Ir);
  W(ib).erg_rho = S.grid(dg).ETA_rho(Ir);

  W(ib).xrg_u   = S.grid(dg).XI_u(Iu);
  W(ib).erg_u   = S.grid(dg).ETA_u(Iu);

  W(ib).xrg_v   = S.grid(dg).XI_v(Iv);
  W(ib).erg_v   = S.grid(dg).ETA_v(Iv);
  
  if (spherical),
    W(ib).Xrg_rho = G(dg).lon_rho(Ir);
    W(ib).Yrg_rho = G(dg).lat_rho(Ir);

    W(ib).Xrg_u   = G(dg).lon_u(Iu);
    W(ib).Yrg_u   = G(dg).lat_u(Iu);

    W(ib).Xrg_v   = G(dg).lon_v(Iv);
    W(ib).Yrg_v   = G(dg).lat_v(Iv);
  else
    W(ib).Xrg_rho = G(dg).x_rho(Ir);
    W(ib).Yrg_rho = G(dg).y_rho(Ir);

    W(ib).Xrg_u   = G(dg).x_u(Iu);
    W(ib).Yrg_u   = G(dg).y_u(Iu);

    W(ib).Xrg_v   = G(dg).x_v(Iv);
    W(ib).Yrg_v   = G(dg).y_v(Iv);
  end
  
  [W(ib).Jrg_rho, W(ib).Irg_rho] = meshgrid(0:Mprg-1, Lprg-1:Lprg-2+Ioffr);
  [W(ib).Jrg_u  , W(ib).Irg_u  ] = meshgrid(0:Mprg-1, Lrg   :Lrg -1+Ioffu);
  [W(ib).Jrg_v  , W(ib).Irg_v  ] = meshgrid(1:Mrg   , Lprg-1:Lprg-2+Ioffv);

  W(ib).Idg_rho   = S.grid(dg).I_rho(Ir);
  W(ib).Jdg_rho   = S.grid(dg).J_rho(Ir);

  W(ib).Idg_u     = S.grid(dg).I_u(Iu);
  W(ib).Jdg_u     = S.grid(dg).J_u(Iu);

  W(ib).Idg_v     = S.grid(dg).I_v(Iv);
  W(ib).Jdg_v     = S.grid(dg).J_v(Iv);

  W(ib).angle     = G(dg).angle(Ir);
  W(ib).f         = G(dg).f(Ir);
  W(ib).h         = G(dg).h(Ir);
  W(ib).pm        = G(dg).pm(Ir);
  W(ib).pn        = G(dg).pn(Ir);

  if (G(dg).curvilinear),
    W(ib).dndx    = G(dg).dndx(Ir);
    W(ib).dmde    = G(dg).dmde(Ir);
  else
    W(ib).dndx    = zeros(size(W(ib).h));
    W(ib).dmde    = zeros(size(W(ib).h));
  end

  W(ib).mask_rho  = G(dg).mask_rho(Ir);
  W(ib).mask_u    = G(dg).mask_u  (Iu);
  W(ib).mask_v    = G(dg).mask_v  (Iv);
end

% Southern boundary of the donor grid lies on the receiver grid northern
% boundary.

if (S.contact(cr).boundary(isouth).okay),
  Joffr = 3;
  Joffu = 3;
  Joffv = 3;

  ib = isouth;

  Jr = S.grid(dg).J_rho <= Joffr;  Jr(:,1) = false;
  Ju = S.grid(dg).J_u   <= Joffu;  Ju(:,1) = false;
  Jv = S.grid(dg).J_v   <= Joffv;

  W(ib).xrg_rho = S.grid(dg).XI_rho(Jr);
  W(ib).erg_rho = S.grid(dg).ETA_rho(Jr);

  W(ib).xrg_u   = S.grid(dg).XI_u(Ju);
  W(ib).erg_u   = S.grid(dg).ETA_u(Ju);

  W(ib).xrg_v   = S.grid(dg).XI_v(Jv);
  W(ib).erg_v   = S.grid(dg).ETA_v(Jv);
  
  if (spherical),
    W(ib).Xrg_rho = G(dg).lon_rho(Jr);
    W(ib).Yrg_rho = G(dg).lat_rho(Jr);
  
    W(ib).Xrg_u   = G(dg).lon_u(Ju);
    W(ib).Yrg_u   = G(dg).lat_u(Ju);

    W(ib).Xrg_v   = G(dg).lon_v(Jv);
    W(ib).Yrg_v   = G(dg).lat_v(Jv);
  else
    W(ib).Xrg_rho = G(dg).x_rho(Jr);
    W(ib).Yrg_rho = G(dg).y_rho(Jr);
  
    W(ib).Xrg_u   = G(dg).x_u(Ju);
    W(ib).Yrg_u   = G(dg).y_u(Ju);

    W(ib).Xrg_v   = G(dg).x_v(Jv);
    W(ib).Yrg_v   = G(dg).y_v(Jv);
  end
  
  [W(ib).Jrg_rho, W(ib).Irg_rho] = meshgrid(Mprg-1:Mprg-2+Joffr, 0:Lpdg-1);
  [W(ib).Jrg_u  , W(ib).Irg_u  ] = meshgrid(Mprg-1:Mprg-2+Joffu, 1:Ldg   );
  [W(ib).Jrg_v  , W(ib).Irg_v  ] = meshgrid(Mrg   :Mrg -1+Joffv, 0:Lpdg-1);

  W(ib).Idg_rho   = S.grid(dg).I_rho(Jr);
  W(ib).Jdg_rho   = S.grid(dg).J_rho(Jr);

  W(ib).Idg_u     = S.grid(dg).I_u(Ju);
  W(ib).Jdg_u     = S.grid(dg).J_u(Ju);

  W(ib).Idg_v     = S.grid(dg).I_v(Jv);
  W(ib).Jdg_v     = S.grid(dg).J_v(Jv);

  W(ib).angle     = G(dg).angle(Jr);
  W(ib).f         = G(dg).f(Jr);
  W(ib).h         = G(dg).h(Jr);
  W(ib).pm        = G(dg).pm(Jr);
  W(ib).pn        = G(dg).pn(Jr);

  if (G(dg).curvilinear),
    W(ib).dndx    = G(dg).dndx(Jr);
    W(ib).dmde    = G(dg).dmde(Jr);
  else
    W(ib).dndx    = zeros(size(W(ib).h));
    W(ib).dmde    = zeros(size(W(ib).h));
  end

  W(ib).mask_rho  = G(dg).mask_rho(Jr);
  W(ib).mask_u    = G(dg).mask_u  (Ju);
  W(ib).mask_v    = G(dg).mask_v  (Jv);
end

% Eastern boundary of the donor grid lies on the receiver grid western
% boundary.

if (S.contact(cr).boundary(ieast).okay),
  Ioffr = 4;
  Ioffu = 4;
  Ioffv = 4;

  ib = ieast;

  Ir = S.grid(dg).I_rho >= -1 + Lpdg - Ioffr;  Ir(end,:) = false;
  Iu = S.grid(dg).I_u   >=  1 + Ldg  - Ioffu;
  Iv = S.grid(dg).I_v   >= -1 + Lpdg - Ioffv;  Iv(end,:) = false;

  W(ib).xrg_rho = S.grid(dg).XI_rho(Ir);
  W(ib).erg_rho = S.grid(dg).ETA_rho(Ir);

  W(ib).xrg_u   = S.grid(dg).XI_u(Iu);
  W(ib).erg_u   = S.grid(dg).ETA_u(Iu);

  W(ib).xrg_v   = S.grid(dg).XI_v(Iv);
  W(ib).erg_v   = S.grid(dg).ETA_v(Iv);
  
  if (spherical),
    W(ib).Xrg_rho = G(dg).lon_rho(Ir);
    W(ib).Yrg_rho = G(dg).lat_rho(Ir);
  
    W(ib).Xrg_u   = G(dg).lon_u(Iu);
    W(ib).Yrg_u   = G(dg).lat_u(Iu);

    W(ib).Xrg_v   = G(dg).lon_v(Iv);
    W(ib).Yrg_v   = G(dg).lat_v(Iv);
  else
    W(ib).Xrg_rho = G(dg).x_rho(Ir);
    W(ib).Yrg_rho = G(dg).y_rho(Ir);
  
    W(ib).Xrg_u   = G(dg).x_u(Iu);
    W(ib).Yrg_u   = G(dg).y_u(Iu);

    W(ib).Xrg_v   = G(dg).x_v(Iv);
    W(ib).Yrg_v   = G(dg).y_v(Iv);
  end
  
  [W(ib).Jrg_rho, W(ib).Irg_rho] = meshgrid(0:Mprg-1,  -Ioffr+1:0);
  [W(ib).Jrg_u  , W(ib).Irg_u  ] = meshgrid(0:Mprg-1, 1-Ioffu+1:1);
  [W(ib).Jrg_v  , W(ib).Irg_v  ] = meshgrid(1:Mrg   ,  -Ioffv+1:0);

  W(ib).Idg_rho   = S.grid(dg).I_rho(Ir);
  W(ib).Jdg_rho   = S.grid(dg).J_rho(Ir);

  W(ib).Idg_u     = S.grid(dg).I_u(Iu);
  W(ib).Jdg_u     = S.grid(dg).J_u(Iu);

  W(ib).Idg_v     = S.grid(dg).I_v(Iv);
  W(ib).Jdg_v     = S.grid(dg).J_v(Iv);
  
  W(ib).angle     = G(dg).angle(Ir);
  W(ib).f         = G(dg).f(Ir);
  W(ib).h         = G(dg).h(Ir);
  W(ib).pm        = G(dg).pm(Ir);
  W(ib).pn        = G(dg).pn(Ir);
 
  if (G(dg).curvilinear),
    W(ib).dndx    = G(dg).dndx(Ir);
    W(ib).dmde    = G(dg).dmde(Ir);
  else
    W(ib).dndx    = zeros(size(W(ib).h));
    W(ib).dmde    = zeros(size(W(ib).h));
  end

  W(ib).mask_rho  = G(dg).mask_rho(Ir);
  W(ib).mask_u    = G(dg).mask_u  (Iu);
  W(ib).mask_v    = G(dg).mask_v  (Iv);
end

% Northern boundary of the donor grid lies on the receiver grid southern
% boundary.

if (S.contact(cr).boundary(inorth).okay),
  Joffr = 4;
  Joffu = 4;
  Joffv = 4;

  ib = inorth;

  Jr = S.grid(dg).J_rho >= -1 + Mpdg - Joffr;  Jr(:,end) = false;
  Ju = S.grid(dg).J_u   >= -1 + Mpdg - Joffu;  Ju(:,end) = false;
  Jv = S.grid(dg).J_v   >=  1 + Mdg  - Joffv;
  
  W(ib).xrg_rho = S.grid(dg).XI_rho(Jr);
  W(ib).erg_rho = S.grid(dg).ETA_rho(Jr);

  W(ib).xrg_u   = S.grid(dg).XI_u(Ju);
  W(ib).erg_u   = S.grid(dg).ETA_u(Ju);

  W(ib).xrg_v   = S.grid(dg).XI_v(Jv);
  W(ib).erg_v   = S.grid(dg).ETA_v(Jv);

  if (spherical),
    W(ib).Xrg_rho = G(dg).lon_rho(Jr);
    W(ib).Yrg_rho = G(dg).lat_rho(Jr);
  
    W(ib).Xrg_u   = G(dg).lon_u(Ju);
    W(ib).Yrg_u   = G(dg).lat_u(Ju);

    W(ib).Xrg_v   = G(dg).lon_v(Jv);
    W(ib).Yrg_v   = G(dg).lat_v(Jv);
  else
    W(ib).Xrg_rho = G(dg).x_rho(Jr);
    W(ib).Yrg_rho = G(dg).y_rho(Jr);
  
    W(ib).Xrg_u   = G(dg).x_u(Ju);
    W(ib).Yrg_u   = G(dg).y_u(Ju);

    W(ib).Xrg_v   = G(dg).x_v(Jv);
    W(ib).Yrg_v   = G(dg).y_v(Jv);
  end

  [W(ib).Jrg_rho, W(ib).Irg_rho] = meshgrid( -Joffr+1:0, 0:Lpdg-1);
  [W(ib).Jrg_u  , W(ib).Irg_u  ] = meshgrid( -Joffu+1:0, 1:Ldg   );
  [W(ib).Jrg_v  , W(ib).Irg_v  ] = meshgrid(1-Joffv+1:1, 0:Lpdg-1);

  W(ib).Idg_rho   = S.grid(dg).I_rho(Jr);
  W(ib).Jdg_rho   = S.grid(dg).J_rho(Jr);

  W(ib).Idg_u     = S.grid(dg).I_u(Ju);
  W(ib).Jdg_u     = S.grid(dg).J_u(Ju);

  W(ib).Idg_v     = S.grid(dg).I_v(Jv);
  W(ib).Jdg_v     = S.grid(dg).J_v(Jv);

  W(ib).angle     = G(dg).angle(Jr);
  W(ib).f         = G(dg).f(Jr);
  W(ib).h         = G(dg).h(Jr);
  W(ib).pm        = G(dg).pm(Jr);
  W(ib).pn        = G(dg).pn(Jr);

  
  if (G(dg).curvilinear),
    W(ib).dndx    = G(dg).dndx(Jr);
    W(ib).dmde    = G(dg).dmde(Jr);
  else
    W(ib).dndx    = zeros(size(W(ib).h));
    W(ib).dmde    = zeros(size(W(ib).h));
  end

  W(ib).mask_rho  = G(dg).mask_rho(Jr);
  W(ib).mask_u    = G(dg).mask_u  (Ju);
  W(ib).mask_v    = G(dg).mask_v  (Jv);
end

%--------------------------------------------------------------------------
%  Convert contact data to vectors.  Currently, it is assumed that the
%  contact region is adjacent to a single boundary. In the future, it
%  is possible to have two or three boundary connections.  If this is
%  case, the contact information is the local work structure, W.  We
%  just need to concatenate all the contact points.
%--------------------------------------------------------------------------

for ib=1:4,
  if (S.contact(cr).boundary(ib).okay),
    C.xrg_rho  = W(ib).xrg_rho;
    C.erg_rho  = W(ib).erg_rho;
    C.xrg_u    = W(ib).xrg_u;
    C.erg_u    = W(ib).erg_u;
    C.xrg_v    = W(ib).xrg_v;
    C.erg_v    = W(ib).erg_v;
    
    C.Xrg_rho  = W(ib).Xrg_rho;
    C.Yrg_rho  = W(ib).Yrg_rho;
    C.Irg_rho  = W(ib).Irg_rho(:);
    C.Jrg_rho  = W(ib).Jrg_rho(:);
    C.Idg_rho  = W(ib).Idg_rho;
    C.Jdg_rho  = W(ib).Jdg_rho;

    C.Xrg_u    = W(ib).Xrg_u;
    C.Yrg_u    = W(ib).Yrg_u;
    C.Irg_u    = W(ib).Irg_u(:);
    C.Jrg_u    = W(ib).Jrg_u(:);
    C.Idg_u    = W(ib).Idg_u;
    C.Jdg_u    = W(ib).Jdg_u;

    C.Xrg_v    = W(ib).Xrg_v;
    C.Yrg_v    = W(ib).Yrg_v;
    C.Irg_v    = W(ib).Irg_v(:);
    C.Jrg_v    = W(ib).Jrg_v(:);
    C.Idg_v    = W(ib).Idg_v;
    C.Jdg_v    = W(ib).Jdg_v;

    C.angle    = W(ib).angle;
    C.f        = W(ib).f;
    C.h        = W(ib).h;
    C.pm       = W(ib).pm;
    C.pn       = W(ib).pn;
    C.dndx     = W(ib).dndx;
    C.dmde     = W(ib).dmde;
    C.mask_rho = W(ib).mask_rho;
    C.mask_u   = W(ib).mask_u;
    C.mask_v   = W(ib).mask_v;

% Remove land contact points.

    if (Lmask),
      indr = C.mask_rho < 0.5;
      indu = C.mask_u   < 0.5;
      indv = C.mask_v   < 0.5;

      C.xrg_rho(indr) = [];
      C.erg_rho(indr) = [];
      C.Xrg_rho(indr) = [];
      C.Yrg_rho(indr) = [];
      C.Irg_rho(indr) = [];
      C.Jrg_rho(indr) = [];
      C.Idg_rho(indr) = [];
      C.Jdg_rho(indr) = [];

      C.xrg_u(indr) = [];
      C.erg_u(indr) = [];
      C.Xrg_u(indu) = [];
      C.Yrg_u(indu) = [];
      C.Irg_u(indu) = [];
      C.Jrg_u(indu) = [];
      C.Idg_u(indu) = [];
      C.Jdg_u(indu) = [];

      C.xrg_v(indr) = [];
      C.erg_v(indr) = [];
      C.Xrg_v(indv) = [];
      C.Yrg_v(indv) = [];
      C.Irg_v(indv) = [];
      C.Jrg_v(indv) = [];
      C.Idg_v(indv) = [];
      C.Jdg_v(indv) = [];

      C.angle   (indr) = [];
      C.f       (indr) = [];
      C.h       (indr) = [];
      C.pm      (indr) = [];
      C.pn      (indr) = [];
      C.dndx    (indr) = [];
      C.dmde    (indr) = [];
      C.mask_rho(indr) = [];
      C.mask_u  (indu) = [];
      C.mask_v  (indv) = [];
    end
  end
end

return

%--------------------------------------------------------------------------

function [C, R] = refinement(cr, dg, rg, Lmask, G, S, MaskInterp)

%
% [P, R] = refinement(cr, dg, rg, Lmask, G, S)
%
% This function sets the contact points in the receiver grid contact
% region using the donor grid data for refinement grids.
%
% On Input:
%
%    cr          Contact region number
%    dg          Donor grid number
%    rg          Receiver grid number
%    Lmask       Switch to remove contact points on land
%    G           Nested grids structure (1 x Ngrid struct array)
%    S           Nested grids information structure (struct array)
%    MaskInterp  Switch to interpolate PSI-, U- and V-mask (true) or
%                  computed from interpolated RHO-mask (false)
%
% On Output:
%
%    C           Contact points structure (struct array):
%
%                  C.xrg_rho(:)            - Receiver grid contact points
%                  C.erg_rho(:)              (XI,ETA) coordinates, (X,Y)
%                  C.Xrg_rho(:)              physical coordinates, and
%                  C.Yrg_rho(:)              and (I,J) indices (RHO-points)
%                  C.Irg_rho(:)              where data is needed from
%                  C.Jrg_rho(:)              donor
%
%                  C.Idg_rho(:)            - Donor (I,J) cell containing
%                  C.Jdg_rho(:)              receiver grid contact point
%
%                  C.xrg_u(:)              - Receiver grid contact points
%                  C.erg_u(:)                (XI,ETA) coordinates, (X,Y)
%                  C.Xrg_u(:)                physical coordinates, and
%                  C.Yrg_u(:)                and (I,J) indices (U-points)
%                  C.Irg_u(:)                where data is needed from
%                  C.Jrg_u(:)                donor
%
%                  C.Idg_u(:)              - Donor (I,J) cell containing
%                  C.Jdg_u(:)                receiver grid contact point
%
%                  C.xrg_v(:)              - Receiver grid contact points
%                  C.erg_v(:)                (XI,ETA) coordinates, (X,Y)
%                  C.Xrg_v(:)                physical coordinates, and
%                  C.Yrg_v(:)                and (I,J) indices (V-points)
%                  C.Irg_v(:)                where data is needed from
%                  C.Jrg_v(:)                donor
%
%                  C.Idg_v(:)              - Donor (I,J) cell containing
%                  C.Jdg_v(:)                receiver grid contact point
%
%                  C.boundary_rho          - Contact point on RHO-boundary
%                  C.boundary_u            - Contact point on   U-boundary
%                  C.boundary_v            - Contact point on   V-boundary
%                                            (initialized to empty)
%
%                                          - Donor data at contact point:
%                  C.angle(:)                  angle between XI-axis & East
%                  C.f(:)                      Coriolis parameter
%                  C.h(:)                      bathymetry
%                  C.pm(:)                     curvilinear  XI-metric, 1/dx
%                  C.pn(:)                     curvilinear ETA-metric, 1/dy
%                  C.dndx(:)                   d(pn)/d(xi)  metric
%                  C.dmde(:)                   d(pm)/d(eta) metric
%                  C.mask_rho(:)               Land/Sea mask at RHO-points
%                  C.mask_u(:)                 Land/Sea mask at U-points
%                  C.mask_v(:)                 land/Sea mask at V-contact
%

% Initialize.

Ngrids      = length(G);             % number of nested grids
Ncontact    = (Ngrids-1)*2;          % number of contact regions
  
method      = 'natural';
spherical   = S.spherical;
Debugging   = true;                  % plot contact points for debugging

C = struct('xrg_rho'     , [], 'erg_rho'     , [],                      ...
           'Xrg_rho'     , [], 'Yrg_rho'     , [],                      ...
           'Irg_rho'     , [], 'Jrg_rho'     , [],                      ...
           'Idg_rho'     , [], 'Jdg_rho'     , [],                      ...
           'xrg_u'       , [], 'erg_u'       , [],                      ...
           'Xrg_u'       , [], 'Yrg_u'       , [],                      ...
           'Irg_u'       , [], 'Jrg_u'       , [],                      ...
           'Idg_u'       , [], 'Jdg_u'       , [],                      ...
           'xrg_v'       , [], 'erg_v'       , [],                      ...
           'Xrg_v'       , [], 'Yrg_v'       , [],                      ...
           'Irg_v'       , [], 'Jrg_v'       , [],                      ...
           'Idg_v'       , [], 'Jdg_v'       , [],                      ...
           'boundary_rho', [], 'boundary_u'  , [], 'boundary_v'  , [],  ...
           'angle'       , [], 'f'           , [], 'h'           , [],  ...
           'pm'          , [], 'pn'          , [],                      ...
           'dndx'        , [], 'dmde'        , [],                      ...
           'mask_rho'    , [], 'mask_u'      , [], 'mask_v'      , []);

% Compute mean grid cell area for donor and receiver grids. In refinement
% the donor and receiver grids have different mean grid cell area.

AreaAvg_dg = mean(mean((1./G(dg).pm) .* (1./G(dg).pn)));
AreaAvg_rg = mean(mean((1./G(rg).pm) .* (1./G(rg).pn)));

%--------------------------------------------------------------------------
% If receiver is a refinement grid (smaller cell area), determine the
% contact points from coarser donor grid (larger cell area).
%--------------------------------------------------------------------------

if ((S.grid(rg).refine_factor > 0) && (AreaAvg_dg > AreaAvg_rg)),

% Extract larger receiver (fine) grid from donor (coarse) grid.

  R = refine_coordinates (cr, dg, rg, G, S, MaskInterp);

% Set the donor grid (Io,Jo) origin coordinates (left-bottom corner) at
% PSI-points used to extract receiver (fine) grid.

  refine_factor = double(S.grid(rg).refine_factor);
  
  delta = 1.0/refine_factor;
  half  = refine_factor/2;
  Io    = min(S.contact(cr).corners.Idg);
  Jo    = min(S.contact(cr).corners.Jdg);

% Set receiver grid contact points: perimeter and outside contact region
% on donor grid.

  if (spherical),

% RHO-contact points.

    [INr,ONr] = inpolygon(R.lon_rho(:), R.lat_rho(:),                   ...
                          S.grid(rg).perimeter.X_psi,                   ...
                          S.grid(rg).perimeter.Y_psi);

    if (any(~INr)),
      C.xrg_rho  = R.xi_rho (~INr);
      C.erg_rho  = R.eta_rho(~INr);
      C.Xrg_rho  = R.lon_rho(~INr);
      C.Yrg_rho  = R.lat_rho(~INr);
      C.Irg_rho  = R.Irg_rho(~INr);
      C.Jrg_rho  = R.Jrg_rho(~INr);
      C.Idg_rho  = floor(Io+(C.Irg_rho-half).*delta);
      C.Jdg_rho  = floor(Jo+(C.Jrg_rho-half).*delta);

      C.angle    = R.angle(~INr);
      C.f        = R.f(~INr);
      C.h        = R.h(~INr);
      C.pm       = R.pm(~INr);
      C.pn       = R.pn(~INr);

      if (G(rg).curvilinear),
        C.dndx   = R.dndx(~INr);
        C.dmde   = R.dmde(~INr);
      else
        C.dndx   = zeros(size(C.h));
        C.dmde   = zeros(size(C.h));
      end

      C.mask_rho = R.mask_rho(~INr);
    end

    if (Debugging),
      figure;
      plot(C.xrg_rho, C.erg_rho, 'r+',                                  ...
           S.grid(dg).XI_rho(:), S.grid(dg).ETA_rho(:), 'co',           ...
           C.Idg_rho+0.5, C.Jdg_rho+0.5, 'bo');
      title(['(XI,ETA) Coordinates, RHO-points: ', blanks(4),           ...
             'Refinement Contact Points (blue circles: Idg,Jdg)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: RHO-points, cyan circles: coarse RHO-points)']);
      ylabel('Boundary Coarse to Fine Cells');
    end

% U-contact points.
    
    [INu,ONu] = inpolygon(R.lon_u(:), R.lat_u(:),                       ...
                          S.grid(rg).perimeter.X_uv,                    ...
                          S.grid(rg).perimeter.Y_uv);

    INu(ONu) = false;                           % add U-points on perimeter

    if (any(~INu)),
      C.xrg_u  = R.xi_u (~INu);
      C.erg_u  = R.eta_u(~INu);
      C.Xrg_u  = R.lon_u(~INu);
      C.Yrg_u  = R.lat_u(~INu);
      C.Irg_u  = R.Irg_u(~INu);
      C.Jrg_u  = R.Jrg_u(~INu);
      C.Idg_u  = ceil(Io+(C.Irg_u-refine_factor).*delta);
      C.Jdg_u  = floor(Jo+(C.Jrg_u-half).*delta);

      C.mask_u = R.mask_u(~INu);
    end

    if (Debugging),
      figure;
      plot(C.xrg_u, C.erg_u, 'r+',                                      ...
           S.grid(dg).XI_u(:), S.grid(dg).ETA_u(:), 'co',               ...
           C.Idg_u, C.Jdg_u+0.5, 'bo');
      title(['(XI,ETA) Coordinates, U-points: ', blanks(4),             ...
             'Refinement Contact Points (blue circles: Idg,Jdg)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: U-points, cyan circles: coarse U-points)']);
      ylabel('Boundary Coarse to Fine Cells');   
    end

% V-contact points.

    [INv,ONv] = inpolygon(R.lon_v(:), R.lat_v(:),                       ...
                          S.grid(rg).perimeter.X_uv,                    ...
                          S.grid(rg).perimeter.Y_uv);

    INv(ONv) = false;                           % add V-points on perimeter
    
    if (any(~INv)),
      C.xrg_v  = R.xi_v (~INv);
      C.erg_v  = R.eta_v(~INv);
      C.Xrg_v  = R.lon_v(~INv);
      C.Yrg_v  = R.lat_v(~INv);
      C.Irg_v  = R.Irg_v(~INv);
      C.Jrg_v  = R.Jrg_v(~INv);
      C.Idg_v  = floor(Io+(C.Irg_v-half).*delta);
      C.Jdg_v  = ceil(Jo+(C.Jrg_v-refine_factor).*delta);

      C.mask_v = R.mask_v(~INv);
    end
  
    if (Debugging),
      figure;
      plot(C.xrg_v, C.erg_v, 'r+',                                      ...
           S.grid(dg).XI_v(:), S.grid(dg).ETA_v(:), 'co',               ...
           C.Idg_v+0.5, C.Jdg_v, 'bo');
      title(['(XI,ETA) Coordinates, V-points: ', blanks(4),             ...
             'Refinement Contact Points (blue circles: Idg,Jdg)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: V-points, cyan circles: coarse V-points)']);
      ylabel('Boundary Coarse to Fine Cells');    
    end

  else

% RHO-contact points.

    [INr,ONr] = inpolygon(R.x_rho(:), R.y_rho(:),                       ...
                          S.grid(rg).perimeter.X_psi,                   ...
                          S.grid(rg).perimeter.Y_psi);

    if (any(~INr)),
      C.xrg_rho  = R.xi_rho (~INr);
      C.erg_rho  = R.eta_rho(~INr);
      C.Xrg_rho  = R.x_rho  (~INr);
      C.Yrg_rho  = R.y_rho  (~INr);
      C.Irg_rho  = R.Irg_rho(~INr);
      C.Jrg_rho  = R.Jrg_rho(~INr);
      C.Idg_rho  = floor(Io+(C.Irg_rho-half).*delta);
      C.Jdg_rho  = floor(Jo+(C.Jrg_rho-half).*delta);

      C.angle    = R.angle(~INr);
      C.f        = R.f(~INr);
      C.h        = R.h(~INr);
      C.pm       = R.pm(~INr);
      C.pn       = R.pn(~INr);

      if (G(rg).curvilinear),
        C.dndx   = R.dndx(~INr);
        C.dmde   = R.dmde(~INr);
      else
        C.dndx   = zeros(size(C.h));
        C.dmde   = zeros(size(C.h));
      end

      C.mask_rho = R.mask_rho(~INr);
    end

    if (Debugging),
      figure;
      plot(C.xrg_rho, C.erg_rho, 'r+',                                  ...
           S.grid(dg).XI_rho(:), S.grid(dg).ETA_rho(:), 'co',           ...
           C.Idg_rho+0.5, C.Jdg_rho+0.5, 'bo');
      title(['(XI,ETA) Coordinates, RHO-points: ', blanks(4),           ...
             'Refinement Contact Points (blue circles: Idg,Jdg)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: RHO-points, cyan circles: coarse RHO-points)']);
      ylabel('Boundary Coarse to Fine Cells');    
    end

% U-contact points.

    [INu,ONu] = inpolygon(R.x_u(:), R.y_u(:),                           ...
                          S.grid(rg).perimeter.X_uv,                    ...
                          S.grid(rg).perimeter.Y_uv);

    INu(ONu) = false;                           % add U-points on perimeter
   
    if (any(~INu)),
      C.xrg_u  = R.xi_u (~INu);
      C.erg_u  = R.eta_u(~INu);
      C.Xrg_u  = R.x_u  (~INu);
      C.Yrg_u  = R.y_u  (~INu);
      C.Irg_u  = R.Irg_u(~INu);
      C.Jrg_u  = R.Jrg_u(~INu);
      C.Idg_u  = ceil(Io+(C.Irg_u-refine_factor).*delta);
      C.Jdg_u  = floor(Jo+(C.Jrg_u-half).*delta);

      C.mask_u = R.mask_u(~INu);
    end

    if (Debugging),
      figure;
      plot(C.xrg_u, C.erg_u, 'r+',                                      ...
           S.grid(dg).XI_u(:), S.grid(dg).ETA_u(:), 'co',               ...
           C.Idg_u, C.Jdg_u+0.5, 'bo');
      title(['(XI,ETA) Coordinates, U-points: ', blanks(4),             ...
             'Refinement Contact Points (blue circles: Idg,Jdg)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: U-points, cyan circles: coarse U-points)']);
      ylabel('Boundary Coarse to Fine Cells');    
    end
    
% V-contact points.
    
    [INv,ONv] = inpolygon(R.x_v(:), R.y_v(:),                           ...
                          S.grid(rg).perimeter.X_uv,                    ...
                          S.grid(rg).perimeter.Y_uv);

    INv(ONv) = false;                           % add V-points on perimeter

    if (any(~INv)),
      C.xrg_v  = R.xi_v (~INv);
      C.erg_v  = R.eta_v(~INv);
      C.Xrg_v  = R.x_v  (~INv);
      C.Yrg_v  = R.y_v  (~INv);
      C.Irg_v  = R.Irg_v(~INv);
      C.Jrg_v  = R.Jrg_v(~INv);
      C.Idg_v  = floor(Io+(C.Irg_v-half).*delta);
      C.Jdg_v  = ceil(Jo+(C.Jrg_v-refine_factor).*delta);

      C.mask_v = R.mask_v(~INv);
    end

    if (Debugging),
      figure;
      plot(C.xrg_v, C.erg_v, 'r+',                                      ...
           S.grid(dg).XI_v(:), S.grid(dg).ETA_v(:), 'co',               ...
           C.Idg_v+0.5, C.Jdg_v, 'bo');
      title(['(XI,ETA) Coordinates, V-points: ', blanks(4),             ...
             'Refinement Contact Points (blue circles: Idg,Jdg)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: V-points, cyan circles: coarse V-points)']);
      ylabel('Boundary Coarse to Fine Cells');    
    end
    
  end

% Remove point on land mask.

  if (Lmask),
    indr = C.mask_rho < 0.5;
    indu = C.mask_u   < 0.5;
    indv = C.mask_v   < 0.5;

    C.xrg_rho(indr) = [];
    C.erg_rho(indr) = [];
    C.Xrg_rho(indr) = [];
    C.Yrg_rho(indr) = [];
    C.Irg_rho(indr) = [];
    C.Jrg_rho(indr) = [];
    C.Idg_rho(indr) = [];
    C.Jdg_rho(indr) = [];

    C.xrg_u(indu) = [];
    C.erg_u(indu) = [];
    C.Xrg_u(indu) = [];
    C.Yrg_u(indu) = [];
    C.Irg_u(indu) = [];
    C.Jrg_u(indu) = [];
    C.Idg_u(indu) = [];
    C.Jdg_u(indu) = [];

    C.xrg_v(indv) = [];
    C.erg_v(indv) = [];
    C.Xrg_v(indv) = [];
    C.Yrg_v(indv) = [];
    C.Irg_v(indv) = [];
    C.Jrg_v(indv) = [];
    C.Idg_v(indv) = [];
    C.Jdg_v(indv) = [];

    C.angle   (indr) = [];
    C.f       (indr) = [];
    C.h       (indr) = [];
    C.pm      (indr) = [];
    C.pn      (indr) = [];
    C.dndx    (indr) = [];
    C.dmde    (indr) = [];
    C.mask_rho(indr) = [];
    C.mask_u  (indu) = [];
    C.mask_v  (indv) = [];
  end
end

%--------------------------------------------------------------------------
% Otherwise if the donor is the finer grid (smaller cell area) and receiver
% the coarser grid (larger cell area), determine coarse grid contact point
% inside the finer grid.  This information will be used for two-way nesting
% for the fine to coarse processing.
%--------------------------------------------------------------------------

if (dg > rg || AreaAvg_rg > AreaAvg_dg),

% Set coaser grid (Io,Jo) origin coordinates (left-bottom corner) at
% PSI-points used to extract finer grid.  Set finer grid center indices
% offset with respect the coarser grid.  Since the corner values may be
% empty or wrong for this contact region, we need to get the (Io,Jo)
% from the conjugate contact region (donor is coarser grid and receiver
% is the finer grid).

  Io     = [];
  Jo     = [];
  delta  = S.grid(dg).refine_factor;
  half   = floor((delta-1)/2);
  offset = double(half+1);

  for my_cr = 1:Ncontact,
    if (S.contact(my_cr).donor_grid    == rg &&                         ...
        S.contact(my_cr).receiver_grid == dg &&                         ...
        S.contact(my_cr).corners.okay),
      Io = min(S.contact(my_cr).corners.Idg);
      Jo = min(S.contact(my_cr).corners.Jdg);
    end
  end

  if (isempty(Io) || isempty(Jo)),
    error(' Unable to determine coaser origin coordinates (Io,Jo)');
  end

% Set receiver (coarse) grid contact points inside perimeter of donor
% (fine) grid. This will be used in the fine two coarse two-way nesting.

  if (spherical),

% RHO-contact points.

    [INr,ONr] = inpolygon(G(rg).lon_rho(:), G(rg).lat_rho(:),           ...
                          S.grid(dg).perimeter.X_psi,                   ...
                          S.grid(dg).perimeter.Y_psi);

    if (any(INr)),
      C.xrg_rho  = S.grid(rg).XI_rho(INr);
      C.erg_rho  = S.grid(rg).ETA_rho(INr);
      C.Xrg_rho  = G(rg).lon_rho(INr);
      C.Yrg_rho  = G(rg).lat_rho(INr);
      C.Irg_rho  = S.grid(rg).I_rho(INr);
      C.Jrg_rho  = S.grid(rg).J_rho(INr);
      C.Idg_rho  = offset+(C.Irg_rho-Io)*delta;
      C.Jdg_rho  = offset+(C.Jrg_rho-Jo)*delta;

      C.angle    = G(rg).angle(INr);
      C.f        = G(rg).f(INr);
      C.h        = G(rg).h(INr);
      C.pm       = G(rg).pm(INr);
      C.pn       = G(rg).pn(INr);

      if (G(rg).curvilinear),
        C.dndx   = G(rg).dndx(INr);
        C.dmde   = G(rg).dmde(INr);
      else
        C.dndx   = zeros(size(C.h));
        C.dmde   = zeros(size(C.h));
      end
      
      C.mask_rho = G(rg).mask_rho(INr);
    end

    if (Debugging),
      figure;
      plot(S.grid(dg).XI_psi,  S.grid(dg).ETA_psi,  'k:',               ...
           S.grid(dg).XI_psi', S.grid(dg).ETA_psi', 'k:',               ...
           S.grid(dg).XI_psi (:,1:delta:end),                           ...
           S.grid(dg).ETA_psi(:,1:delta:end),  'k-',                    ...
           S.grid(dg).XI_psi (1:delta:end,:)',                          ...
           S.grid(dg).ETA_psi(1:delta:end,:)', 'k-',                    ...
           S.grid(dg).XI_rho(:), S.grid(dg).ETA_rho(:), 'r+',           ...
           C.Idg_rho+0.5, C.Jdg_rho+0.5,'bo');
      title(['(XI,ETA) Coordinates, RHO-points: ', blanks(4),           ...
             'two-way Coarse Grid Points (blue circles)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: RHO-points, mesh: PSI-points)']);
      ylabel('Interior Fine To Coarse Cells');
    end

% U-contact points.

    [INu,ONu] = inpolygon(G(rg).lon_u(:), G(rg).lat_u(:),               ...
                          S.grid(dg).perimeter.X_uv,                    ...
                          S.grid(dg).perimeter.Y_uv);

    INu(ONu) = false;                   % remove U-points on perimeter

    if (any(INu)),
      C.xrg_u  = S.grid(rg).XI_u(INu);
      C.erg_u  = S.grid(rg).ETA_u(INu);
      C.Xrg_u  = G(rg).lon_u(INu);
      C.Yrg_u  = G(rg).lat_u(INu);
      C.Irg_u  = S.grid(rg).I_u(INu);
      C.Jrg_u  = S.grid(rg).J_u(INu);
      C.Idg_u  = 1+(C.Irg_u-Io)*delta;
      C.Jdg_u  = offset+(C.Jrg_u-Jo)*delta;

      C.mask_u = G(rg).mask_u(INu);
    end

    if (Debugging),
      figure;
      plot(S.grid(dg).XI_v,  S.grid(dg).ETA_v,  'k:',                   ...
           S.grid(dg).XI_v', S.grid(dg).ETA_v', 'k:',                   ...
           S.grid(dg).XI_v (offset+1:end-offset,1:delta:end),           ...
           S.grid(dg).ETA_v(offset+1:end-offset,1:delta:end), 'k-',     ...
           S.grid(dg).XI_v (offset+1:delta:end,:)',                     ...
           S.grid(dg).ETA_v(offset+1:delta:end,:)', 'k-',               ...
           S.grid(dg).perimeter.XI_psi,                                 ...
           S.grid(dg).perimeter.ETA_psi, 'ks',                          ...
           S.grid(dg).XI_u(:), S.grid(dg).ETA_u(:), 'r+',               ...
           C.Idg_u, C.Jdg_u+0.5,'bo');
      title(['(XI,ETA) Coordinates, U-points: ', blanks(4),             ...
             'two-way Coarse Grid Points (blue circles)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: U-points, black squares: PSI-points, ',       ...
              'mesh: V-points)']);
      ylabel('Interior Fine To Coarse Cells');
    end

% V-contact points.

    [INv,ONv] = inpolygon(G(rg).lon_v(:), G(rg).lat_v,                  ...
                          S.grid(dg).perimeter.X_uv,                    ...
                          S.grid(dg).perimeter.Y_uv);

    INv(ONv) = false;                   % remove V-points on perimeter

    if (any(INv)),
      C.xrg_v  = S.grid(rg).XI_v(INv);
      C.erg_v  = S.grid(rg).ETA_v(INv);
      C.Xrg_v  = G(rg).lon_v(INv);
      C.Yrg_v  = G(rg).lat_v(INv);
      C.Irg_v  = S.grid(rg).I_v(INv);
      C.Jrg_v  = S.grid(rg).J_v(INv);
      C.Idg_v  = offset+(C.Irg_v-Io).*delta;
      C.Jdg_v  = 1+(C.Jrg_v-Jo).*delta;
      
      C.mask_v = G(rg).mask_v(INv);
    end

    if (Debugging),
      figure;
      plot(S.grid(dg).XI_u,  S.grid(dg).ETA_u,  'k:',                   ...
           S.grid(dg).XI_u', S.grid(dg).ETA_u', 'k:',                   ...
           S.grid(dg).XI_u (:,offset+1:delta:end),                      ...
           S.grid(dg).ETA_u(:,offset+1:delta:end), 'k-',                ...
           S.grid(dg).XI_u (1:delta:end,offset+1:end-offset)',          ...
           S.grid(dg).ETA_u(1:delta:end,offset+1:end-offset)', 'k-',    ...
           S.grid(dg).perimeter.XI_psi,                                 ...
           S.grid(dg).perimeter.ETA_psi, 'ks',                          ...
           S.grid(dg).XI_v(:), S.grid(dg).ETA_v(:), 'r+',               ...
           C.Idg_v+0.5, C.Jdg_v,'bo');
      title(['(XI,ETA) Coordinates, V-points: ', blanks(4),             ...
             'two-way Coarse Grid Points (blue circles)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: V-points, black squares: PSI-points, ',       ...
              'mesh: U-points)']);
      ylabel('Interior Fine To Coarse Cells');
    end
  
  else

% RHO-contact points.

    [INr,ONr] = inpolygon(G(rg).x_rho(:), G(rg).y_rho(:),               ...
                          S.grid(dg).perimeter.X_psi,                   ...
                          S.grid(dg).perimeter.Y_psi);

    if (any(INr)),
      C.xrg_rho  = S.grid(rg).XI_rho(INr);
      C.erg_rho  = S.grid(rg).ETA_rho(INr);
      C.Xrg_rho  = G(rg).x_rho(INr);
      C.Yrg_rho  = G(rg).y_rho(INr);
      C.Irg_rho  = S.grid(rg).I_rho(INr);
      C.Jrg_rho  = S.grid(rg).J_rho(INr);
      C.Idg_rho  = offset+(C.Irg_rho-Io)*delta;
      C.Jdg_rho  = offset+(C.Jrg_rho-Jo)*delta;
    
      C.angle    = G(rg).angle(INr);
      C.f        = G(rg).f(INr);
      C.h        = G(rg).h(INr);
      C.pm       = G(rg).pm(INr);
      C.pn       = G(rg).pn(INr);

      if (G(rg).curvilinear),
        C.dndx   = G(rg).dndx(INr);
        C.dmde   = G(rg).dmde(INr);
      else
        C.dndx   = zeros(size(C.h));
        C.dmde   = zeros(size(C.h));
      end

      C.mask_rho = G(rg).mask_rho(INr);
    end

    if (Debugging),
      figure;
      plot(S.grid(dg).XI_psi,  S.grid(dg).ETA_psi,  'k:',               ...
           S.grid(dg).XI_psi', S.grid(dg).ETA_psi', 'k:',               ...
           S.grid(dg).XI_psi (:,1:delta:end),                           ...
           S.grid(dg).ETA_psi(:,1:delta:end),  'k-',                    ...
           S.grid(dg).XI_psi (1:delta:end,:)',                          ...
           S.grid(dg).ETA_psi(1:delta:end,:)', 'k-',                    ...
           S.grid(dg).XI_rho(:), S.grid(dg).ETA_rho(:), 'r+',           ...
           C.Idg_rho+0.5, C.Jdg_rho+0.5,'bo');
      title(['(XI,ETA) Coordinates, RHO-points: ', blanks(4),           ...
             'two-way Coarse Grid Points (blue circles)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: RHO-points, mesh: PSI-points)']);
      ylabel('Interior Fine To Coarse Cells');
    end
    
% U-contact points.

    [INu,ONu] = inpolygon(G(rg).x_u(:), G(rg).y_u(:),                   ...
                          S.grid(dg).perimeter.X_uv,                    ...
                          S.grid(dg).perimeter.Y_uv);

    INu(ONu) = false;                   % remove U-points on perimeter

    if (any(INu)),
      C.xrg_u  = S.grid(rg).XI_u(INu);
      C.erg_u  = S.grid(rg).ETA_u(INu);
      C.Xrg_u  = G(rg).x_u(INu);
      C.Yrg_u  = G(rg).y_u(INu);
      C.Irg_u  = S.grid(rg).I_u(INu);
      C.Jrg_u  = S.grid(rg).J_u(INu);
      C.Idg_u  = 1+(C.Irg_u-Io)*delta;
      C.Jdg_u  = offset+(C.Jrg_u-Jo)*delta;

      C.mask_u = G(rg).mask_u(INu);
    end

    if (Debugging),
      figure;
      plot(S.grid(dg).XI_v,  S.grid(dg).ETA_v,  'k:',                   ...
           S.grid(dg).XI_v', S.grid(dg).ETA_v', 'k:',                   ...
           S.grid(dg).XI_v (offset+1:end-offset,1:delta:end),           ...
           S.grid(dg).ETA_v(offset+1:end-offset,1:delta:end), 'k-',     ...
           S.grid(dg).XI_v (offset+1:delta:end,:)',                     ...
           S.grid(dg).ETA_v(offset+1:delta:end,:)', 'k-',               ...
           S.grid(dg).perimeter.XI_psi,                                 ...
           S.grid(dg).perimeter.ETA_psi, 'ks',                          ...
           S.grid(dg).XI_u(:), S.grid(dg).ETA_u(:), 'r+',               ...
           C.Idg_u, C.Jdg_u+0.5,'bo');
      title(['(XI,ETA) Coordinates, U-points: ', blanks(4),             ...
             'two-way Coarse Grid Points (blue circles)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: U-points, black squares: PSI-points, ',       ...
              'mesh: V-points)']);
      ylabel('Interior Fine To Coarse Cells');
    end
    
% V-contact points.

    [INv,ONv] = inpolygon(G(rg).x_v(:), G(rg).y_v(:),                   ...
                          S.grid(dg).perimeter.X_uv,                    ...
                          S.grid(dg).perimeter.Y_uv);

    INv(ONv) = false;                   % remove V-points on perimeter

    if (any(INv)),
      C.xrg_v  = S.grid(rg).XI_v(INv);
      C.erg_v  = S.grid(rg).ETA_v(INv);
      C.Xrg_v  = G(rg).x_v(INv);
      C.Yrg_v  = G(rg).y_v(INv);
      C.Irg_v  = S.grid(rg).I_v(INv);
      C.Jrg_v  = S.grid(rg).J_v(INv);
      C.Idg_v  = offset+(C.Irg_v-Io).*delta;
      C.Jdg_v  = 1+(C.Jrg_v-Jo).*delta;

      C.mask_v = G(rg).mask_v(INv);
    end
    
    if (Debugging),
      figure;
      plot(S.grid(dg).XI_u,  S.grid(dg).ETA_u,  'k:',                   ...
           S.grid(dg).XI_u', S.grid(dg).ETA_u', 'k:',                   ...
           S.grid(dg).XI_u (:,offset+1:delta:end),                      ...
           S.grid(dg).ETA_u(:,offset+1:delta:end), 'k-',                ...
           S.grid(dg).XI_u (1:delta:end,offset+1:end-offset)',          ...
           S.grid(dg).ETA_u(1:delta:end,offset+1:end-offset)', 'k-',    ...
           S.grid(dg).perimeter.XI_psi,                                 ...
           S.grid(dg).perimeter.ETA_psi, 'ks',                          ...
           S.grid(dg).XI_v(:), S.grid(dg).ETA_v(:), 'r+',               ...
           C.Idg_v+0.5, C.Jdg_v,'bo');
      title(['(XI,ETA) Coordinates,  V-points: ', blanks(4),            ...
             'two-way Coarse Grid Points (blue circles)']);
      xlabel(['cr = ', num2str(cr), blanks(4),                          ...
              'dg = ', num2str(dg), blanks(4),                          ...
              'rg = ', num2str(rg), blanks(8),                          ...
              '(red plus: V-points, black squares: PSI-points, ',       ...
              'mesh: U-points)']);
      ylabel('Interior Fine To Coarse Cells');
    end
    
  end

% Set intermediate refinement structure to empty.

  R.spherical = [];
  R.uniform   = [];

  R.xi_rho  = [];
  R.eta_rho = [];

  R.x_rho   = [];
  R.y_rho   = [];
  if (spherical),
    R.lon_rho = [];
    R.lat_rho = [];
  end

  R.Irg_rho = [];
  R.Jrg_rho = [];

  R.xi_psi  = [];
  R.eta_psi = [];

  R.x_psi   = [];
  R.y_psi   = [];
  if (spherical),
    R.lon_psi = [];
    R.lat_psi = [];
  end

  R.Irg_psi = [];
  R.Jrg_psi = [];

  R.xi_u  = [];
  R.eta_u = [];

  R.x_u   = [];
  R.y_u   = [];
  if (spherical),
    R.lon_u = [];
    R.lat_u = [];
  end

  R.Irg_u = [];
  R.Jrg_u = [];

  R.xi_v  = [];
  R.eta_v = [];
  
  R.x_v   = [];
  R.y_v   = [];
  if (spherical),
    R.lon_v = [];
    R.lat_v = [];
  end

  R.Irg_v = [];
  R.Jrg_v = [];

  R.angle    = [];
  R.f        = [];
  R.h        = [];
  R.pm       = [];
  R.pn       = [];
  R.dndx     = [];
  R.dmde     = [];
  R.mask_rho = [];
  R.mask_psi = [];
  R.mask_u   = [];
  R.mask_v   = [];

end

return

%--------------------------------------------------------------------------

function S = boundary_contact(Sinp)

%
% S = boundary_contact(Sinp)
%
% This function determines which contact points lay on the reciever grid
% boundary. This information is used in ROMS to avoid applying lateral
% boundary conditions at those points.
%
% On Input:
%
%    Sinp       Nested grids information structure (struct array)
%
% On Output:
%
%    S          Sinp structure plus lateral boundary switch (struct array)
%
% This function adds the following field to the S structure for each
% contact region, cr:
%
%    S.contact(cr).point.boundary_rho      - Contact point on RHO-boundary
%    S.contact(cr).point.boundary_u        - Contact point on   U-boundary
%    S.contact(cr).point.boundary_v        - Contact point on   V-boundary
%

% Initialize output structure with input structure.
  
S = Sinp;

UsePolygon = true;                   %  Use "inpolygon" fiunction

%--------------------------------------------------------------------------
% Determine which contact point lay on receiver grid boundary.
%--------------------------------------------------------------------------

for cr=1:S.Ncontact,

  rg = S.contact(cr).receiver_grid;
  Lm = S.grid(rg).L-1;
  Mm = S.grid(rg).M-1;
  
% Contact points on RHO-boundary.

  if (UsePolygon)
    [IN,ON] = inpolygon(S.contact(cr).point.Xrg_rho,                    ...
                        S.contact(cr).point.Yrg_rho,                    ...
                        S.grid(rg).perimeter.X_psi,                     ...
                        S.grid(rg).perimeter.Y_psi);
  else
    NC = length(S.contact(cr).point.Xrg_rho);

    ON = false(size(S.contact(cr).point.Xrg_rho));

    Xper = S.grid(rg).perimeter.X_psi;
    Yper = S.grid(rg).perimeter.Y_psi;
  
    for nc=1:NC,
      Xrg_rho = S.contact(cr).point.Xrg_rho(nc);
      Yrg_rho = S.contact(cr).point.Yrg_rho(nc);
    
      ind = find(abs(Xper-Xrg_rho) < 4*eps &                            ...
                 abs(Yper-Yrg_rho) < 4*eps);
      if (~isempty(ind)),
        ON(nc) = true;
      end
    end
  end
    
  S.contact(cr).point.boundary_rho = ON;

% Contact points on U-boundary.  

  if (UsePolygon)
    [IN,ON] = inpolygon(S.contact(cr).point.Xrg_u,                      ...
                        S.contact(cr).point.Yrg_u,                      ...
                        S.grid(rg).perimeter.X_uv,                      ...
                        S.grid(rg).perimeter.Y_uv);
  else
    NC = length(S.contact(cr).point.Xrg_u);

    ON = false(size(S.contact(cr).point.Xrg_u));

    Xper = S.grid(rg).perimeter.X_uv;
    Yper = S.grid(rg).perimeter.Y_uv;
  
    for nc=1:NC,
      Xrg_u = S.contact(cr).point.Xrg_u(nc);
      Yrg_u = S.contact(cr).point.Yrg_u(nc);
    
      ind = find(abs(Xper-Xrg_u) < 4*eps &                              ...
                 abs(Yper-Yrg_u) < 4*eps);
      if (~isempty(ind)),
        ON(nc) = true;
      end
    end
  end

  S.contact(cr).point.boundary_u = ON;

% Contact points on V-boundary.

  if (UsePolygon)
    [IN,ON] = inpolygon(S.contact(cr).point.Xrg_v,                      ...
                        S.contact(cr).point.Yrg_v,                      ...
                        S.grid(rg).perimeter.X_uv,                      ...
                        S.grid(rg).perimeter.Y_uv);
  else
    NC = length(S.contact(cr).point.Xrg_v);

    ON = false(size(S.contact(cr).point.Xrg_v));

    Xper = S.grid(rg).perimeter.X_uv;
    Yper = S.grid(rg).perimeter.Y_uv;
  
    for nc=1:NC,
      Xrg_v = S.contact(cr).point.Xrg_v(nc);
      Yrg_v = S.contact(cr).point.Yrg_v(nc);
    
      ind = find(abs(Xper-Xrg_v) < 4*eps &                              ...
                 abs(Yper-Yrg_v) < 4*eps);
      if (~isempty(ind)),
        ON(nc) = true;
      end
    end
  end

  S.contact(cr).point.boundary_v = ON;

end

return

%--------------------------------------------------------------------------

function S = Hweights(G, Sinp, ImposeMask)

%
% S = Hweights(G, Sinp)
%
% This function computes the contact points horizontal interpolation
% weights to set data in the receiver grid contact region using donor
% grid data.
%
% On Input:
%
%    G          Nested grids fields structure (struct array)
%
%    Sinp       Nested grids information structure (struct array)
%
%    ImposeMask Switch to scale weights with the land/sea mask.
%
% On Output:
%
%    S          Sinp structure plus horizontal interpolation weights
%                 (struct array)
%
% This function add the following fields to the S structure for each
% contact region, cr:
%
%    S.Lweights(cr).H_rho(4,:)            - Contact points linear weights
%    S.Lweights(cr).H_u  (4,:)              (H) to horizontally interpolate
%    S.Lweights(cr).H_v  (4,:)              reciever data from donor grid
%     
%    S.Qweights(cr).H_rho(9,:)            - Contact points quadratic weights
%    S.Qweights(cr).H_u  (9,:)              (H) to horizontally interpolate
%    S.Qweights(cr).H_v  (9,:)              reciever data from donor grid

% Initialize output structure with input structure.
  
S = Sinp;

% Initialize

Ncontact  = S.Ncontact;
spherical = S.spherical;

Ldebug = true;

%--------------------------------------------------------------------------
% Set horizontal interpolation weights.
%--------------------------------------------------------------------------

for cr=1:Ncontact,

  dg = S.contact(cr).donor_grid;
  rg = S.contact(cr).receiver_grid;

% RHO-contact points.  Recall that we need to shift I and J by one since
% Matlat does not support zero-indices.

  if (S.contact(cr).refinement && (dg > rg)),
    Hzero = zeros(size(S.contact(cr).point.Irg_rho));

    S.Lweights(cr).H_rho(1,:) = Hzero;
    S.Lweights(cr).H_rho(2,:) = Hzero;
    S.Lweights(cr).H_rho(3,:) = Hzero;
    S.Lweights(cr).H_rho(4,:) = Hzero;

    S.Qweights(cr).H_rho(1,:) = Hzero;
    S.Qweights(cr).H_rho(2,:) = Hzero;
    S.Qweights(cr).H_rho(3,:) = Hzero;
    S.Qweights(cr).H_rho(4,:) = Hzero;
    S.Qweights(cr).H_rho(5,:) = Hzero;
    S.Qweights(cr).H_rho(6,:) = Hzero;
    S.Qweights(cr).H_rho(7,:) = Hzero;
    S.Qweights(cr).H_rho(8,:) = Hzero;
    S.Qweights(cr).H_rho(9,:) = Hzero;
    
    Rcompute = false;
  else
    [Ir,Jr] = size(S.grid(dg).I_rho);
    RindexB = sub2ind([Ir, Jr],                                         ...
                      S.contact(cr).point.Idg_rho+1,                    ...
                      min(Jr,S.contact(cr).point.Jdg_rho));  % (Idg, Jdg-1)
    RindexO = sub2ind([Ir, Jr],                                         ...
                      S.contact(cr).point.Idg_rho+1,                    ...
                      S.contact(cr).point.Jdg_rho+1);        % (Idg, Jdg)
    RindexT = sub2ind([Ir, Jr],                                         ...
                      S.contact(cr).point.Idg_rho+1,                    ...
                      min(Jr,S.contact(cr).point.Jdg_rho+2));% (Idg, Jdg+1)
    RindexR = sub2ind([Ir, Jr],                                         ...
                      min(Ir, S.contact(cr).point.Idg_rho+2),           ...
                      S.contact(cr).point.Jdg_rho+1);        % (Idg+1, Jdg)
    Rcompute = true;
  end
  if (Rcompute),
    if (S.contact(cr).refinement)
      pr = (S.contact(cr).point.xrg_rho - S.grid(dg).XI_rho(RindexO))./ ...
           (S.grid(dg).XI_rho(RindexR)  - S.grid(dg).XI_rho(RindexO));
      qr = (S.contact(cr).point.erg_rho - S.grid(dg).ETA_rho(RindexO))./...
           (S.grid(dg).ETA_rho(RindexT) - S.grid(dg).ETA_rho(RindexO));
    else
      if (spherical),
        pr = (S.contact(cr).point.Xrg_rho - G(dg).lon_rho(RindexO)) ./  ...
             (G(dg).lon_rho(RindexR) - G(dg).lon_rho(RindexO));
        qr = (S.contact(cr).point.Yrg_rho - G(dg).lat_rho(RindexO)) ./  ...
             (G(dg).lat_rho(RindexT) - G(dg).lat_rho(RindexO));
      else
        pr = (S.contact(cr).point.Xrg_rho - G(dg).x_rho(RindexO)) ./    ...
             (G(dg).x_rho(RindexR) - G(dg).x_rho(RindexO));
        qr = (S.contact(cr).point.Yrg_rho - G(dg).y_rho(RindexO)) ./    ...
             (G(dg).y_rho(RindexT) - G(dg).y_rho(RindexO));
      end
    end

    S.Lweights(cr).H_rho = linear_weights(pr, qr,                       ...
                                          RindexO, RindexR, RindexT,    ...
                                          ImposeMask, G(dg).mask_rho);

    S.Qweights(cr).H_rho = quadratic_weights(1-pr, qr,                  ...
                                             RindexB, RindexO, RindexT, ...
                                             S.grid(rg).refine_factor,  ...
                                             ImposeMask, G(dg).mask_rho);

  end

% U-contact points.  Recall that we need to shift I by one since
% Matlat does not support zero-indices.

  if (S.contact(cr).refinement && (dg > rg)),
    Hzero = zeros(size(S.contact(cr).point.Irg_u));

    S.Lweights(cr).H_u(1,:) = Hzero;
    S.Lweights(cr).H_u(2,:) = Hzero;
    S.Lweights(cr).H_u(3,:) = Hzero;
    S.Lweights(cr).H_u(4,:) = Hzero;

    S.Qweights(cr).H_u(1,:) = Hzero;
    S.Qweights(cr).H_u(2,:) = Hzero;
    S.Qweights(cr).H_u(3,:) = Hzero;
    S.Qweights(cr).H_u(4,:) = Hzero;
    S.Qweights(cr).H_u(5,:) = Hzero;
    S.Qweights(cr).H_u(6,:) = Hzero;
    S.Qweights(cr).H_u(7,:) = Hzero;
    S.Qweights(cr).H_u(8,:) = Hzero;
    S.Qweights(cr).H_u(9,:) = Hzero;
    
    Ucompute = false;
  else
    [Iu,Ju] = size(S.grid(dg).I_u);
    UindexB = sub2ind([Iu, Ju],                                         ...
                      S.contact(cr).point.Idg_u  ,                      ...
                      min(Ju, S.contact(cr).point.Jdg_u));   % (Idg, Jdg-1)
    UindexO = sub2ind([Iu, Ju],                                         ...
                      S.contact(cr).point.Idg_u  ,                      ...
                      min(Ju, S.contact(cr).point.Jdg_u+1)); % (Idg, Jdg)
    UindexT = sub2ind([Iu, Ju],                                         ...
                      S.contact(cr).point.Idg_u  ,                      ...
                      min(Ju, S.contact(cr).point.Jdg_u+2)); % (Idg, Jdg+1)
    UindexR = sub2ind([Iu, Ju],                                         ...
                      min(Iu, S.contact(cr).point.Idg_u+1),             ...
                      min(Ju, S.contact(cr).point.Jdg_u+1)); % (Idg+1, Jdg)
    Ucompute = true;
  end
  if (Ucompute)    
    if (S.contact(cr).refinement)
      pu = (S.contact(cr).point.xrg_u - S.grid(dg).XI_u(UindexO))./     ...
           (S.grid(dg).XI_u(UindexR)  - S.grid(dg).XI_u(UindexO));
      qu = (S.contact(cr).point.erg_u - S.grid(dg).ETA_u(UindexO))./    ...
           (S.grid(dg).ETA_u(UindexT) - S.grid(dg).ETA_u(UindexO));
    else
      if (spherical),
        pu = (S.contact(cr).point.Xrg_u - G(dg).lon_u(UindexO)) ./      ...
             (G(dg).lon_u(UindexR) - G(dg).lon_u(UindexO));
        qu = (S.contact(cr).point.Yrg_u - G(dg).lat_u(UindexO)) ./      ...
             (G(dg).lat_u(UindexT) - G(dg).lat_u(UindexO)) ;
      else
        pu = (S.contact(cr).point.Xrg_u - G(dg).x_u(UindexO)) ./        ...
             (G(dg).x_u(UindexR) - G(dg).x_u(UindexO));
        qu = (S.contact(cr).point.Yrg_u - G(dg).y_u(UindexO)) ./        ...
             (G(dg).y_u(UindexT) - G(dg).y_u(UindexO));
      end
    end
      
    S.Lweights(cr).H_u = linear_weights(pu, qu,                         ...
                                        UindexO, UindexR, UindexT,      ...
                                        ImposeMask, G(dg).mask_u);

    S.Qweights(cr).H_u = quadratic_weights(1-pu, qu,                    ...
                                           UindexB, UindexO, UindexT,   ...
                                           S.grid(rg).refine_factor,    ...
                                           ImposeMask, G(dg).mask_u);

  end

% V-contact points.  Recall that we need to shift J by one since
% Matlat does not support zero-indices.

  if (S.contact(cr).refinement && (dg > rg)),
    Hzero = zeros(size(S.contact(cr).point.Irg_v));
    S.Lweights(cr).H_v(1,:) = Hzero;
    S.Lweights(cr).H_v(2,:) = Hzero;
    S.Lweights(cr).H_v(3,:) = Hzero;
    S.Lweights(cr).H_v(4,:) = Hzero;

    S.Qweights(cr).H_v(1,:) = Hzero;
    S.Qweights(cr).H_v(2,:) = Hzero;
    S.Qweights(cr).H_v(3,:) = Hzero;
    S.Qweights(cr).H_v(4,:) = Hzero;
    S.Qweights(cr).H_v(5,:) = Hzero;
    S.Qweights(cr).H_v(6,:) = Hzero;
    S.Qweights(cr).H_v(7,:) = Hzero;
    S.Qweights(cr).H_v(8,:) = Hzero;
    S.Qweights(cr).H_v(9,:) = Hzero;

    Vcompute = false;
  else
    [Iv,Jv] = size(S.grid(dg).I_v);
    VindexB = sub2ind([Iv, Jv],                                         ...
                      min(Iv, S.contact(cr).point.Idg_v+1),             ...
                      min(Jv, S.contact(cr).point.Jdg_v-1)); % (Idg, Jdg-1)
    VindexO = sub2ind([Iv, Jv],                                         ...
                      min(Iv, S.contact(cr).point.Idg_v+1),             ...
                      S.contact(cr).point.Jdg_v  );          % (Idg, Jdg)
    VindexT = sub2ind([Iv, Jv],                                         ...
                      min(Iv, S.contact(cr).point.Idg_v+1),             ...
                      min(Jv, S.contact(cr).point.Jdg_v+1)); % (Idg, Jdg+1)
    VindexR = sub2ind([Iv, Jv],                                         ...
                      min(Iv, S.contact(cr).point.Idg_v+2),             ...
                      S.contact(cr).point.Jdg_v  );          % (Idg+1, Jdg)
    Vcompute = true;
  end
  if (Vcompute),
    if (S.contact(cr).refinement)
      pv = (S.contact(cr).point.xrg_v - S.grid(dg).XI_v(VindexO))./     ...
           (S.grid(dg).XI_v(VindexR)  - S.grid(dg).XI_v(VindexO));
      qv = (S.contact(cr).point.erg_v - S.grid(dg).ETA_v(VindexO))./    ...
           (S.grid(dg).ETA_v(VindexT) - S.grid(dg).ETA_v(VindexO));
    else
      if (spherical),
        pv = (S.contact(cr).point.Xrg_v - G(dg).lon_v(VindexO)) ./      ...
             (G(dg).lon_v(VindexR) - G(dg).lon_v(VindexO));
        qv = (S.contact(cr).point.Yrg_v - G(dg).lat_v(VindexO)) ./      ...
             (G(dg).lat_v(VindexT) - G(dg).lat_v(VindexO)) ;
      else
        pv = (S.contact(cr).point.Xrg_v - G(dg).x_v(VindexO)) ./        ...
             (G(dg).x_v(VindexR) - G(dg).x_v(VindexO));
        qv = (S.contact(cr).point.Yrg_v - G(dg).y_v(VindexO)) ./        ...
             (G(dg).y_v(VindexT) - G(dg).y_v(VindexO));
      end
    end
      
    S.Lweights(cr).H_v = linear_weights(pv, qv,                         ...
                                        VindexO, VindexR, VindexT,      ...
                                        ImposeMask, G(dg).mask_v);

    S.Qweights(cr).H_v = quadratic_weights(1-pv, qv,                    ...
                                           VindexB, VindexO, VindexT,   ...
                                           S.grid(rg).refine_factor,    ...
                                           ImposeMask, G(dg).mask_v);
  end

% Debugging.

  if (Ldebug),

    frmt = ['%7.7i %12.5e %12.5e %12.5e %12.5e %12.5e %5i %7i %7i ',    ...
              '%10.5f %10.5f %10.5f %5i %7i %10.5f %10.5f %10.5f\n'];
    if (S.contact(cr).refinement),
      Labels = ['   n     Hweight(1)   Hweight(2)   Hweight(3) ',       ...
                '  Hweight(4)       SUM       Idg   indxO   indxR  ',   ...
                '  XIdgO       XIrg      XIdgR     Jdg   indxT  ',      ...
                '  ETAdgO     ETArg      ETAdgT'];
    else
      if (spherical)
        Labels = ['   n     Hweight(1)   Hweight(2)   Hweight(3) ',     ...
                  '  Hweight(4)       SUM       Idg   indxO ',          ...
                  '  indxR    dgLonO      Xrg       dgLonR  ',          ...
                  '  Jdg   indxT    dgLatO      Yrg       dgLatT'];
      else
        Labels = ['   n     Hweight(1)   Hweight(2)   Hweight(3) ',     ...
                  '  Hweight(4)       SUM       Idg   indxO ',          ...
                  '  indxR     XdgO       Xrg        XdgR   ',          ...
                  '  Jdg   indxT     YdgO       Yrg        YdxT'];
      end
    end
    
    if (cr == 1),
      Rname = 'Rweigths.txt';
      Rout = fopen(Rname, 'w+');
      if (Rout < 0),
        error(['Hweights - Cannot create debugging file: ', Rname, '.']);
      end
      s = 'Interpolation Weights at RHO-contact points: ';
      fprintf (Rout, '%s donor = %2.2i, receiver = %2.2i \n\n', s, dg, rg);
    else
      fprintf (Rout, ' \n\n');
      fprintf (Rout, 'Contact Region = %2.2i\n\n', cr);
    end
    
    if (Rcompute),
      fprintf (Rout, '%s \n\n', Labels);
      if (S.contact(cr).refinement),
        for n=1:length(S.contact(cr).point.Xrg_rho),
          fprintf (Rout, frmt,                                          ...
                   [n,                                                  ...
                    transpose(S.Lweights(cr).H_rho(:,n)),               ...
                    sum(S.Lweights(cr).H_rho(:,n)),                     ...
                    S.contact(cr).point.Idg_rho(n)+1,                   ...
                    RindexO(n),                                         ...
                    RindexR(n),                                         ...
                    S.grid(dg).XI_rho(RindexO(n)),                      ...
                    S.contact(cr).point.xrg_rho(n),                     ...
                    S.grid(dg).XI_rho(RindexR(n)),                      ...
                    S.contact(cr).point.Jdg_rho(n)+1,                   ...
                    RindexT(n),                                         ...
                    S.grid(dg).ETA_rho(RindexO(n)),                     ...
                    S.contact(cr).point.erg_rho(n),                     ...
                    S.grid(dg).ETA_rho(RindexT(n))]);
        end
      else
        if (spherical)
          for n=1:length(S.contact(cr).point.Xrg_rho),
            fprintf (Rout, frmt,                                        ...
                     [n,                                                ...
                      transpose(S.Lweights(cr).H_rho(:,n)),             ...
                      sum(S.Lweights(cr).H_rho(:,n)),                   ...
                      S.contact(cr).point.Idg_rho(n)+1,                 ...
                      RindexO(n),                                       ...
                      RindexR(n),                                       ...
                      G(dg).lon_rho(RindexO(n)),                        ...
                      S.contact(cr).point.Xrg_rho(n),                   ...
                      G(dg).lon_rho(RindexR(n)),                        ...
                      S.contact(cr).point.Jdg_rho(n)+1,                 ...
                      RindexT(n),                                       ...
                      G(dg).lat_rho(RindexO(n)),                        ...
                      S.contact(cr).point.Yrg_rho(n),                   ...
                      G(dg).lat_rho(RindexT(n))]);
          end
        else
          for n=1:length(S.contact(cr).point.Xrg_rho),
            fprintf (Rout, frmt,                                        ...
                     [n,                                                ...
                      transpose(S.Lweights(cr).H_rho(:,n)),             ...
                      sum(S.Lweights(cr).H_rho(:,n)),                   ...
                      S.contact(cr).point.Idg_rho(n)+1,                 ...
                      RindexO(n),                                       ...
                      RindexR(n),                                       ...
                      G(dg).x_rho(RindexO(n)),                          ...
                      S.contact(cr).point.Xrg_rho(n),                   ...
                      G(dg).x_rho(RindexR(n)),                          ...
                      S.contact(cr).point.Jdg_rho(n)+1,                 ...
                      RindexT(n),                                       ...
                      G(dg).y_rho(RindexO(n)),                          ...
                      S.contact(cr).point.Yrg_rho(n),                   ...
                      G(dg).y_rho(RindexT(n))]);
          end
        end
      end
    else
      fprintf (Rout,'  Weights are zero because they are not needed.\n\n');
    end

    if (cr == 1),
      Uname = 'Uweigths.txt';
      Uout = fopen(Uname, 'w+');
      if (Uout < 0),
        error(['Hweights - Cannot create debugging file: ', Uname, '.']);
      end
      s = 'Interpolation Weights at U-contact points: ';
      fprintf (Uout, '%s donor = %2.2i, receiver = %2.2i \n\n', s, dg, rg);
    else
      fprintf (Uout, ' \n\n');
      fprintf (Uout, 'Contact Region = %2.2i\n\n', cr);
    end

    if (Ucompute),
      fprintf (Uout, '%s \n\n', Labels);
      if (S.contact(cr).refinement),
        for n=1:length(S.contact(cr).point.Xrg_u),
          fprintf (Uout, frmt,                                          ...
                   [n,                                                  ...
                    transpose(S.Lweights(cr).H_u(:,n)),                 ...
                    sum(S.Lweights(cr).H_u(:,n)),                       ...
                    S.contact(cr).point.Idg_u(n)+1,                     ...
                    UindexO(n),                                         ...
                    UindexR(n),                                         ...
                    S.grid(dg).XI_u(UindexO(n)),                        ...
                    S.contact(cr).point.xrg_u(n),                       ...
                    S.grid(dg).XI_u(UindexR(n)),                        ...
                    S.contact(cr).point.Jdg_u(n)+1,                     ...
                    UindexT(n),                                         ...
                    S.grid(dg).ETA_u(UindexO(n)),                       ...
                    S.contact(cr).point.erg_u(n),                       ...
                    S.grid(dg).ETA_u(UindexT(n))]);
        end
      else
        if (spherical)
          for n=1:length(S.contact(cr).point.Xrg_u),
            fprintf (Uout, frmt,                                        ...
                     [n,                                                ...
                      transpose(S.Lweights(cr).H_u(:,n)),               ...
                      sum(S.Lweights(cr).H_u(:,n)),                     ...
                      S.contact(cr).point.Idg_u(n)+1,                   ...
                      UindexO(n),                                       ...
                      UindexR(n),                                       ...
                      G(dg).lon_u(UindexO(n)),                          ...
                      S.contact(cr).point.Xrg_u(n),                     ...
                      G(dg).lon_u(UindexR(n)),                          ...
                      S.contact(cr).point.Jdg_u(n)+1,                   ...
                      UindexT(n),                                       ...
                      G(dg).lat_u(UindexO(n)),                          ...
                      S.contact(cr).point.Yrg_u(n),                     ...
                      G(dg).lat_u(UindexT(n))]);
          end
        else
          for n=1:length(S.contact(cr).point.Xrg_u),
            fprintf (Uout, frmt,                                        ...
                     [n,                                                ...
                      transpose(S.Lweights(cr).H_u(:,n)),               ...
                      sum(S.Lweights(cr).H_u(:,n)),                     ...
                      S.contact(cr).point.Idg_u(n)+1,                   ...
                      UindexO(n),                                       ...
                      UindexR(n),                                       ...
                      G(dg).x_u(UindexO(n)),                            ...
                      S.contact(cr).point.Xrg_u(n),                     ...
                      G(dg).x_u(UindexR(n)),                            ...
                      S.contact(cr).point.Jdg_u(n)+1,                   ...
                      UindexT(n),                                       ...
                      G(dg).y_u(UindexO(n)),                            ...
                      S.contact(cr).point.Yrg_u(n),                     ...
                      G(dg).y_u(UindexT(n))]);
          end
        end
      end
    else
      fprintf (Uout,'  Weights are zero because they are not needed.\n\n');
    end

    if (cr == 1),
      Vname = 'Vweigths.txt';
      Vout = fopen(Vname, 'w+');
      if (Vout < 0),
        error(['Hweights - Cannot create debugging file: ', Vname, '.']);
      end
      s = 'Interpolation Weights at V-contact points: ';
      fprintf (Vout, '%s donor = %2.2i, receiver = %2.2i \n\n', s, dg, rg);
    else
      fprintf (Vout, ' \n\n');
      fprintf (Vout, 'Contact Region = %2.2i\n\n', cr);    
    end

    if (Vcompute),
      fprintf (Vout, '%s \n\n', Labels);
      if (S.contact(cr).refinement),
        for n=1:length(S.contact(cr).point.Xrg_v),
          fprintf (Vout, frmt,                                          ...
                   [n,                                                  ...
                    transpose(S.Lweights(cr).H_v(:,n)),                 ...
                    sum(S.Lweights(cr).H_v(:,n)),                       ...
                    S.contact(cr).point.Idg_v(n)+1,                     ...
                    VindexO(n),                                         ...
                    VindexR(n),                                         ...
                    S.grid(dg).XI_v(VindexO(n)),                        ...
                    S.contact(cr).point.xrg_v(n),                       ...
                    S.grid(dg).XI_v(VindexR(n)),                        ...
                    S.contact(cr).point.Jdg_v(n)+1,                     ...
                    VindexT(n),                                         ...
                    S.grid(dg).ETA_v(VindexO(n)),                       ...
                    S.contact(cr).point.erg_v(n),                       ...
                    S.grid(dg).ETA_v(VindexT(n))]);
        end
      else
        if (spherical)
          for n=1:length(S.contact(cr).point.Xrg_v),
            fprintf (Vout, frmt,                                        ...
                     [n,                                                ...
                      transpose(S.Lweights(cr).H_v(:,n)),               ...
                      sum(S.Lweights(cr).H_v(:,n)),                     ...
                      S.contact(cr).point.Idg_v(n)+1,                   ...
                      VindexO(n),                                       ...
                      VindexR(n),                                       ...
                      G(dg).lon_v(VindexO(n)),                          ...
                      S.contact(cr).point.Xrg_v(n),                     ...
                      G(dg).lon_v(VindexR(n)),                          ...
                      S.contact(cr).point.Jdg_v(n)+1,                   ...
                      VindexT(n),                                       ...
                      G(dg).lat_v(VindexO(n)),                          ...
                      S.contact(cr).point.Yrg_v(n),                     ...
                      G(dg).lat_v(VindexT(n))]);
          end
        else
          for n=1:length(S.contact(cr).point.Xrg_v),
            fprintf (Vout, frmt,                                        ...
                     [n,                                                ...
                      transpose(S.Lweights(cr).H_v(:,n)),               ...
                      sum(S.Lweights(cr).H_v(:,n)),                     ...
                      S.contact(cr).point.Idg_v(n)+1,                   ...
                      VindexO(n),                                       ...
                      VindexR(n),                                       ...
                      G(dg).x_v(VindexO(n)),                            ...
                      S.contact(cr).point.Xrg_v(n),                     ...
                      G(dg).x_v(VindexR(n)),                            ...
                      S.contact(cr).point.Jdg_v(n)+1,                   ...
                      VindexT(n),                                       ...
                      G(dg).y_v(VindexO(n)),                            ...
                      S.contact(cr).point.Yrg_v(n),                     ...
                      G(dg).y_v(VindexT(n))]);
          end
        end
      end
    else
      fprintf (Vout,'  Weights are zero because they are not needed.\n\n');
    end
  end
  clear RindexO RindexR RindexT Rcompute
  clear UindexO UindexR UindexT Ucompute
  clear VindexO VindexR VindexT Vcompute
end

if (Ldebug),
  fclose (Rout);
  fclose (Uout);
  fclose (Vout);
end  

return

%--------------------------------------------------------------------------

function W = linear_weights(p, q, index1, index2, index4, ImposeMask, mask)

%
% W = linear_weights(p, q, index1, index2, index4, ImposeMask, mask)
%
% This function computes the contact points horizontal linear
% interpolation weights given the fractional distances (p,q) from
% the donor grid.  The weights are adjusted in the presence of
% land/sea masking.
%
% On Input:
%
%    p          Contact points fractional I-distance with respect
%                 the donor grid cell (vector)
%
%    q          Contact points fractional J-distance with respect
%                 the donor grid cell (vector)
%
%    index1     Single linear indexes for 2D subscripts in donor
%                 grid cell corner 1 (Idg,Jdg) for each contact
%                 point (vector).
%
%    index2     Single linear indexes for 2D subscripts in donor
%                 grid cell corner 2 (Idg+1,Jdg) for each contact
%                 point (vector).
%
%    index4     Single linear indexes for 2D subscripts in donor
%                 grid cell corner 4 (Idg,Jdg+1) for each contact
%                 point (vector).
%
%    ImposeMask Switch to scale weights with the land/sea mask.
%
%    mask       Donor grid land/sea masking (2D array).
%
% On Output:
%
%    W(1:4,:)   Contact point linear interpolation weights (2D array):
%
%                 W(1,:)    interpolation weight for corner 1
%                 W(2,:)    interpolation weight for corner 2
%                 W(3,:)    interpolation weight for corner 3
%                 W(4,:)    interpolation weight for corner 4
%
%
% The following diagrams show the conventions for linear interpolation:
%
%                 index4         index3
%
%                   4______________3  (Idg+1,Jdg+1)
%                   |  :           |
%                   |  :       p   |
%                   | 1-q   <----->|
%                   |  :           |
%                 Jr|..:....x   :  |
%                   |       .   :  |
%                   |  1-p  .   :  |
%                   |<------>   q  |
%                   |       .   :  |
%                   |       .   :  |
%                   |___________:__|
%         (Idg,Jdg) 1       Ir     2
%
%                 index1         index2
%
% For linear interpolation and all water points, the weights are:
%
%     W(1,:) = (1 - p) * (1 - q)
%     W(2,:) = p * (1 - q)
%     W(3,:) = p * q
%     W(4,:) = (1 - p) * q
%
% In ROMS, the linear interpolation is carried out as:
%
%     V(Irg, Jrg) = W(1,:) * F2D(Idg,  Jdg  ) +
%                   W(2,:) * F2D(Idg+1,Jdg  ) + 
%                   W(3,:) * F2D(Idg+1,Jdg+1) +
%                   W(4,:) * F2D(Idg,  Jdg+1) 
%

% Initalize.

Npoints = length(p);
index3  = index4 + 1;                % (Idg+1, Jdg+1)

% If any of the input fractional distances is NaN, it implies that
% we were diving by zero in the calling function. Therefore, the
% contact point is concident to the donor grid physical grid boundary
% and there are not values at either Idg+1 or Jdg+1.

if (any(isnan(p))),
  p(isnan(p)) = 0;
end

if (any(isnan(q))),
  q(isnan(qr)) = 0;
end

%--------------------------------------------------------------------------
% Compute linear interpolation weigths.
%--------------------------------------------------------------------------

W  = zeros([4 Npoints]);
LW = [0 0 0 0];

for n=1:Npoints,
  LW(1) = mask(index1(n)) * (1 - p(n)) * (1 - q(n));
  LW(2) = mask(index2(n)) * p(n) * (1 - q(n));
  LW(3) = mask(index3(n)) * p(n) * q(n);
  LW(4) = mask(index4(n)) * (1 -p(n)) * q(n);

  if (ImposeMask),
    MaskSum = mask(index1(n)) + mask(index2(n)) +                       ...
              mask(index3(n)) + mask(index4(n));

    if (MaskSum < 4),                % at least one of the corners is land
      LWsum = sum(LW);
      if (LWsum > 0),
        LW = LW ./ LWsum;            % using only water points
      else
        LW(1:4) = 0;                 % all donor grid corners are on land
      end
    end
  end
    
  W(:,n) = LW;
end

% Impose positive zero.

ind0 = find(abs(W) < 100*eps);
if (~isempty(ind0)),
  W(ind0) = 0;
end

% Impose exact unity.

ind1 = find(abs(W - 1) < 100*eps);
if (~isempty(ind1)),
  W(ind1) = 1;
end

return

%--------------------------------------------------------------------------

function W = quadratic_weights(p, q, index2, index5, index8, rfactor,   ...
                               ImposeMask, mask)

%
% W = quadratic_weights(p, q, index2, index5, index8, rfactor,
%                       ImposeMask, mask)
%
% This function computes the contact points horizontal quadratic
% interpolation weights given the fractional distances (p,q) from
% the donor grid.  The weights are adjusted in the presence of
% land/sea masking.
%
% If refined grids (rfactor > 0), the quadratic interpolation weights are
% conservative.  That is, the coarse-to-fine and fine-to-coarse are
% reversible (Clark and Farley, 1984):
%
%        SUM(F(:,:,...) / rfactor^2 = C(Idg,Jdg,...)
%
% On Input:     (see diagram below)
%
%    p          Contact points fractional I-distance with respect
%                 the donor grid point 5 (vector)
%
%    q          Contact points fractional J-distance with respect
%                 the donor grid point 5 (vector)
%
%    index2     Single linear indexes for 2D subscripts in donor
%                 grid cell corner 2 (Idg,Jdg-1) for each contact
%                 point (vector).
%
%    index5     Single linear indexes for 2D subscripts in donor
%                 grid cell corner 5 (Idg,Jdg) for each contact
%                 point (vector).
%
%    index8     Single linear indexes for 2D subscripts in donor
%                 grid cell corner 8 (Idg,Jdg+1) for each contact
%                 point (vector).
%
%    rfactor    Refine grid factor (scalar)
%
%                 rfactor = 0      quadratic interpolation
%                                  (composite/mosaic grids)
%
%                 rfactor > 0      conservative quadratic interpolation
%                                  (refine grids)
%
%    ImposeMask Switch to scale weights with the land/sea mask.
%
%    mask       Donor grid land/sea masking (2D array).
%
% On Output:
%
%    W(1:9,:)   Contact point quadratic interpolation weights (2D array):
%
%                 W(1,:)    interpolation weight in terms of donor point 1
%                 W(2,:)    interpolation weight in terms of donor point 2
%                 W(3,:)    interpolation weight in terms of donor point 3
%                 W(4,:)    interpolation weight in terms of donor point 4
%                 W(5,:)    interpolation weight in terms of donor point 5
%                 W(6,:)    interpolation weight in terms of donor point 6
%                 W(7,:)    interpolation weight in terms of donor point 7
%                 W(8,:)    interpolation weight in terms of donor point 8
%                 W(9,:)    interpolation weight in terms of donor point 9
%
%
% The following diagrams show the conventions for quadratic interpolation:
%
%               index7           index8           index9
%
%                 Rm               Ro               Rp
%
%            Jd+1 7________________8________________9     Sp
%                 |                |                |
%                 |                |                |
%                 |                |<---p--->       |
%                 |                |                |
%                 |             Jr |........x   :   |
%                 |                |       .    :   |
%                 |                |       .    q   |
%                 |                |       .    :   |
%            Jd   |________________|____________:___|
%                 4                5       Ir       6     So
%                 |                |                |
%               index4           index5           index6
%                 |                |                |
%                 |                |                |
%                 |                |                |
%                 |                |                |
%                 |                |                |
%            Jd-1 |________________|________________|
%                 1                2                3     Sm
%
%                Id-1              Id              Id+1  
%
%               index1           index2           index3
%
% For quadratic interpolation the coefficients are:
%
%     Rm = 0.5 * p * (p - 1) + alpha
%     Ro = (1 - p^2)         - 2 * alpha
%     Rp = 0.5 * p * (p + 1) + alpha
%
%     Sm = 0.5 * q * (q - 1) + alpha
%     So = (1 - q^2)         - 2 * alpha
%     Sp = 0.5 * q * (q + 1) + alpha
%  
% where      alpha = [(1/rfactor)^2 - 1] / 24       for rfactor > 0
%                                                   (conservative)
%
%            alpha = 0                              for rfactor = 0
%
% The quadratic interpolation weights for all water points are:
%
%     W(1,:) = Rm * Sm
%     W(2,:) = Ro * Sm
%     W(3,:) = Rp * Sm
%     W(4,:) = Rm * So
%     W(5,:) = Ro * So
%     W(6,:) = Rp * So
%     W(7,:) = Rm * Sp
%     W(8,:) = Ro * Sp
%     W(9,:) = Rp * Sp
%
% In ROMS, the quadratic interpolation is carried out as:
%
%     F(Irg, Jrg) = W(1,:) * C(Idg-1, Jdg-1) +
%                   W(2,:) * C(Idg  , Jdg-1) +
%                   W(3,:) * C(Idg+1, Jdg-1) +
%                   W(4,:) * C(Idg-1, Jdg  ) +
%                   W(5,:) * C(Idg  , Jdg  ) +
%                   W(6,:) * C(Idg+1, Jdg  ) +
%                   W(7,:) * C(Idg-1, Jdg+1) +
%                   W(8,:) * C(Idg  , Jdg+1) +
%                   W(9,:) * C(Idg+1, Jdg+1)
%

% Initalize.

Npoints = length(p);
index1  = index2 - 1;                % (Idg-1, Jdg-1)
index3  = index2 + 1;                % (Idg+1, Jdg-1)
index4  = index5 - 1;                % (Idg-1, Jdg  )
index6  = index5 + 1;                % (Idg+1, Jdg  )
index7  = index8 - 1;                % (Idg-1, Jdg+1)
index9  = index8 + 1;                % (Idg+1, Jdg+1)

% Compute coefficient for conservative quadratic interpolation.

if (rfactor > 0)
  alpha = ((1 / rfactor)^2 - 1) / 24;
else
  alpha = 0;
end

% If any of the input fractional distances is NaN, it implies that
% we were diving by zero in the calling function. Therefore, the
% contact point is concident to the donor grid physical grid boundary
% and there are not values at either Idg+1 or Jdg+1.

if (any(isnan(p))),
  p(isnan(p)) = 0;
end

if (any(isnan(q))),
  q(isnan(qr)) = 0;
end

%--------------------------------------------------------------------------
% Compute quadratic interpolation weigths.
%--------------------------------------------------------------------------

W  = zeros([9 Npoints]);
QW = zeros([1 9]);

Rm = 0.5 .* p .* (p - 1) + alpha;
Ro = (1 - p .* p)        - 2 * alpha;
Rp = 0.5 .* p .* (p + 1) + alpha;

Sm = 0.5 .* q .* (q - 1) + alpha;
So = (1 - q .* q)        - 2 * alpha;
Sp = 0.5 .* q .* (q + 1) + alpha;

for n=1:Npoints,
  QW(1) = mask(index1(n)) * Rm(n) * Sm(n);
  QW(2) = mask(index2(n)) * Ro(n) * Sm(n);
  QW(3) = mask(index3(n)) * Rp(n) * Sm(n);
  QW(4) = mask(index4(n)) * Rm(n) * So(n);
  QW(5) = mask(index5(n)) * Ro(n) * So(n);
  QW(6) = mask(index6(n)) * Rp(n) * So(n);
  QW(7) = mask(index7(n)) * Rm(n) * Sp(n);
  QW(8) = mask(index8(n)) * Ro(n) * Sp(n);
  QW(9) = mask(index9(n)) * Rp(n) * Sp(n);

  if (ImposeMask),
    MaskSum = mask(index1(n)) + mask(index2(n)) +  mask(index3(n)) +    ...
              mask(index4(n)) + mask(index5(n)) +  mask(index6(n)) +    ...
              mask(index7(n)) + mask(index8(n)) +  mask(index9(n));

    if (MaskSum < 9),          % at least one of the donor points is land
      QWsum = sum(QW);
      if (QWsum > 0),
        QW = QW ./ QWsum;      % using only water points
      else
        QW(1:9) = 0;           % all donor points are on land
      end
    end
  end
  
  W(:,n) = QW;    
end

% Impose positive zero.

ind0 = find(abs(W) < 100*eps);
if (~isempty(ind0)),
  W(ind0) = 0;
end

% Impose exact unity.

ind1 = find(abs(W - 1) < 100*eps);
if (~isempty(ind1)),
  W(ind1) = 1;
end

return