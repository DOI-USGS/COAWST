function S = grid_connections(G, Sinp)

%
% GRID_CONNECTIONS:  Sets Nested Grids Connectivity
%
% S = grid_connection(G, Sinp)
%
% This function appends the nested grid connectivity fields between
% donor and receiver grids for each contact region to the Nested
% Grids Structure, Sinp.
%
% On Input:
%
%    G          Nested Grids Structure (1 x Ngrids struct array)
%
%                 G(ng) = get_roms_grid ( char(Gnames(ng)) )
%
%    Sinp       Contact Points Structure (struct array)
%                 See "contact.m" for full list of fields
%
% On Output:
%
%    S          Updated nested grids contact points structure
%                 (struct array)
%
% This function adds the following fields to the S structure for each
% contact region, cr:
%
%    S.contact(cr).donor_grid              - Donor grid number
%    S.contact(cr).receiver_grid           - Receiver grid number
%    S.contact(cr).coincident              - Coincident boundary switch 
%    S.contact(cr).composite               - Composite grid switch
%    S.contact(cr).hybrid                  - hybrid nested grids switch
%    S.contact(cr).mosaic                  - Mosaic grid switch
%    S.contact(cr).refinement              - Refinement grid switch
%
%    S.contact(cr).interior.okay           - true/false logical
%    S.contact(cr).interior.Xdg(:)         - (X,Y) coordinates and (I,J)
%    S.contact(cr).interior.Ydg(:)           indices of donor grid points
%    S.contact(cr).interior.Idg(:)           inside the receiver grid
%    S.contact(cr).interior.Jdg(:)           perimeter, [] if false 
%
%    S.contact(cr).corners.okay            - true/false logical
%    S.contact(cr).corners.Xdg(:)          - (X,Y) coordinates and (I,J)
%    S.contact(cr).corners.Ydg(:)            indices of donor grid points
%    S.contact(cr).corners.Idg(:)            corners laying on receiver
%    S.contact(cr).corners.Idg(:)            grid perimeter, [] if false
%
%    S.contact(cr).boundary(ib).okay       - true/false logical
%    S.contact(cr).boundary(ib).match(:)   - Donor matching points logical
%    S.contact(cr).boundary(ib).Xdg(:)     - (X,Y) coordinates and (I,J)
%    S.contact(cr).boundary(ib).Ydg(:)       indices of donor boundary
%    S.contact(cr).boundary(ib).Idg(:)       points laying on receiver
%    S.contact(cr).boundary(ib).Jdg(:)       grid perimeter, [] if false
%
% If the receiver is a refinement grid, the "interior" sub-structure
% contains the donor grid RHO-points outside the receiver grid perimeter.
% Otherwise, it contains the donor grid RHO-points inside the receiver
% grid perimenter.  
%

% svn $Id: grid_connections.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize output structure with input structure.
  
S = Sinp;

% Initialize.

iwest  = S.western_edge;         % western  boundary edge index
isouth = S.southern_edge;        % southern boundary edge index
ieast  = S.eastern_edge;         % eastern  boundary edge index
inorth = S.northern_edge;        % northern boundary edge index

B(1:4) = struct('ind', []);      % boundary structure

adjacent  = [ieast, inorth, iwest, isouth];   % Receiver grid boundary
spherical = S.spherical;                      % spherical grid switch

debugging = false;
%debugging = true;  %jcw

% Compute grid maximum spacing.  In refinement, the donor grid is
% coarser than receiver grid so the grid spacing is larger.

DXmax = zeros([1 S.Ngrids]);
DYmax = zeros([1 S.Ngrids]);

for ng=1:S.Ngrids
  DXmax(ng) = max(1./G(ng).pm(:));
  DYmax(ng) = max(1./G(ng).pn(:));
end

% Determine which grids are inside of others but itself.  The structure
% C stores the direct and indirect connections between grids.

C(1:S.Ngrids) = struct('ToGrid', [], 'count', 0);

IsInside = false([S.Ngrids S.Ngrids]);

for dg=1:S.Ngrids
  for rg=1:S.Ngrids

    if (dg ~= rg)

      if (spherical)    
        [IN,~] = inpolygon(G(dg).lon_psi, G(dg).lat_psi,                ...
                           S.grid(rg).perimeter.X_psi,                  ...
                           S.grid(rg).perimeter.Y_psi);    
        if (debugging)
          figure;
          plot(G(dg).lon_psi(:),                                        ...
               G(dg).lat_psi(:),           'b+',                        ...
               S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 'ro');
        end
      else    
        [IN,~] = inpolygon(G(dg).x_psi, G(dg).y_psi,                    ...
                           S.grid(rg).perimeter.X_psi,                  ...
                           S.grid(rg).perimeter.Y_psi);
        if (debugging)
          figure;
          plot(G(dg).x_psi(:),                                          ...
               G(dg).y_psi(:),             'b+',                        ...
               S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 'ro');
        end
      end

      if any(IN)   %orig
%     if any(IN(:))   %jcw
%      if (sum(IN(:))==size(IN,1).*size(IN,2))
        IsInside(rg,dg) = true;
        C(dg).ToGrid = [C(dg).ToGrid rg];      
        C(dg).count  = C(dg).count + 1;      
      end
    
      if (debugging)
        str = {'F', 'T'};
        title(['dg = ', num2str(dg), '  (blue)', blanks(4),             ...
               'rg = ', num2str(rg), '  (red)' , blanks(4),             ...
               'IsInside(rg,dg) = ', char(str(1+IsInside(rg,dg)))]);
      end
    
    end
  end
end

% If telescoping grids, there is more than one connection between grids.
% Reject indirect connections.

for dg = 1:S.Ngrids
  if (C(dg).count > 1)
    rg = max(C(dg).ToGrid);             % direct connection
    for i = 1: C(dg).count
      ng = C(dg).ToGrid(i);
      if (ng ~= rg)
        IsInside(ng,dg) = false;        % reject indirect connection   %jcw commented out
      end
    end
  end
end

% Determine if telescoping grids exist. In ROMS, a telescoping grid
% is a refined grid (refine_factor > 0) containing a finer grid
% inside. Notice that under this definition the coarser grid (ng=1,
% refine_factor=0) is not considered a telescoping grid. The last
% cascading grid is not considered a telescoping grid since it does
% not contains a finer grid inside.  The telescoping definition is
% a technical one, and it facilitates the algorithm design.

CoarserDonor = zeros([1 S.Ngrids]);           % finer grid coarser donor
FinerDonor   = zeros([1 S.Ngrids]);           % coarser grid finer donor
RefinedGrid  = [S.grid.refine_factor] > 0;    % refinement grid switch
Telescoping  = false([1 S.Ngrids]);           % telescoping grid switch

if any(RefinedGrid)
  RefinedGrid(1) = true;
end

for dg=1:S.Ngrids
  for rg=1:S.Ngrids

    if (dg ~= rg)
      if (RefinedGrid(rg) && IsInside(dg,rg) &&                         ...
          S.grid(rg).refine_factor &&                                   ...
          DXmax(dg) > DXmax(rg) && DYmax(dg) > DYmax(rg))
        CoarserDonor(rg) = dg;
        FinerDonor(dg) = rg;
      end

      if (CoarserDonor(rg) == dg && S.grid(dg).refine_factor > 0)
        Telescoping(dg)=true;
      end
    end
  
  end
end

%--------------------------------------------------------------------------
% Loop over all contact regions.
%--------------------------------------------------------------------------

cr = 0;

connected = false([S.Ngrids S.Ngrids]);

disp(blanks(2));
disp(' Summary of Contact Regions Processed:')
disp(blanks(2));
for ng=1:S.Ngrids
  disp(['   Grid ',num2str(ng, '%2.2i'), ': ', G(ng).grid_name]);
end
disp(blanks(2));
disp('   Contact   Donor   Receiver');
disp('    Region    Grid       Grid');
disp(blanks(2));

for dg=1:S.Ngrids
  for rg=1:S.Ngrids
    if (dg ~= rg)

      contact.donor_grid    = dg;             % donor grid number
      contact.receiver_grid = rg;             % receiver grid number
      contact.coincident    = false;
      contact.composite     = false;
      contact.hybrid        = false;
      contact.mosaic        = false;

      if (S.grid(dg).refine_factor > 0 || S.grid(rg).refine_factor > 0)
        contact.refinement = true;
      else
        contact.refinement = false;
      end  

% Determine which points from donor grid are inside the reciever grid.
% If the receiver is a refinment grid, consider only donor RHO-points
% outside the receiver grid. Otherwise, consider only points inside
% the receiver grid.

      if (spherical)
        [IN,~] = inpolygon(G(dg).lon_rho, G(dg).lat_rho,                ...
                           S.grid(rg).perimeter.X_psi,                  ...
                           S.grid(rg).perimeter.Y_psi);
        if (any(IN(:)))
          X = G(dg).lon_rho(:);
          Y = G(dg).lat_rho(:);
          I = S.grid(dg).I_rho(:);
          J = S.grid(dg).J_rho(:);

          contact.interior.okay = true;
          if (S.grid(rg).refine_factor > 0)
            contact.interior.Xdg = X(~IN(:));
            contact.interior.Ydg = Y(~IN(:));
            contact.interior.Idg = I(~IN(:));
            contact.interior.Jdg = J(~IN(:));
          else
            contact.interior.Xdg = X(IN(:));
            contact.interior.Ydg = Y(IN(:));
            contact.interior.Idg = I(IN(:));
            contact.interior.Jdg = J(IN(:));
          end
          clear I J X Y
        else
          contact.interior.okay  = false;
          contact.interior.Xdg   = [];
          contact.interior.Ydg   = [];
          contact.interior.Idg   = [];
          contact.interior.Jdg   = [];
        end
      else
        [IN,~] = inpolygon(G(dg).x_rho, G(dg).y_rho,                    ...
                           S.grid(rg).perimeter.X_psi,                  ...
                           S.grid(rg).perimeter.Y_psi);
        if (any(IN(:)))
          X = G(dg).x_rho(:);
          Y = G(dg).y_rho(:);
          I = S.grid(dg).I_rho(:);
          J = S.grid(dg).J_rho(:);

          contact.interior.okay = true;
          if (S.grid(rg).refine_factor > 0)
            contact.interior.Xdg = X(~IN(:));
            contact.interior.Ydg = Y(~IN(:));
            contact.interior.Idg = I(~IN(:));
            contact.interior.Jdg = J(~IN(:));
          else
            contact.interior.Xdg = X(IN(:));
            contact.interior.Ydg = Y(IN(:));
            contact.interior.Idg = I(IN(:));
            contact.interior.Jdg = J(IN(:));
          end
          clear I J X Y
        else
          contact.interior.okay  = false;
          contact.interior.Xdg   = [];
          contact.interior.Ydg   = [];
          contact.interior.Idg   = [];
          contact.interior.Jdg   = [];
        end
      end

% Determine if any of the donor grid domain corners lay on the receiver
% grid perimeter.
  
      [~,ON] = inpolygon(S.grid(dg).corners.X,                          ...
                         S.grid(dg).corners.Y,                          ...
                         S.grid(rg).perimeter.X_psi,                    ...
                         S.grid(rg).perimeter.Y_psi);
      if (any(ON))
        contact.corners.okay = true;
        myindex = S.grid(dg).corners.index;
        if (spherical)
          contact.corners.Xdg = G(dg).lon_psi(myindex(ON));
          contact.corners.Ydg = G(dg).lat_psi(myindex(ON));
        else
          contact.corners.Xdg = G(dg).x_psi(myindex(ON));
          contact.corners.Ydg = G(dg).y_psi(myindex(ON));
        end
        contact.corners.Idg   = S.grid(dg).I_psi(myindex(ON));
        contact.corners.Jdg   = S.grid(dg).J_psi(myindex(ON));
      else
        contact.corners.okay  = false;
        contact.corners.Xdg   = [];
        contact.corners.Ydg   = [];
        contact.corners.Idg   = [];
        contact.corners.Jdg   = [];
      end

% If receiver (rg) is a refinement grid and donor (dg) is the course
% source grid, determine which donor grid points lay in the receiver
% grid perimeter.

      if ((S.grid(rg).refine_factor > 0) &&                             ...
          (CoarserDonor(rg) == dg || FinerDonor(rg) == dg) &&           ...
          contact.interior.okay)

        set_boundary = true;

        if (spherical)
          [~,ON] = inpolygon(G(dg).lon_psi, G(dg).lat_psi,              ...
                             S.grid(rg).perimeter.X_psi,                ...
                             S.grid(rg).perimeter.Y_psi);

          I = S.grid(dg).I_psi(ON); Is = min(I);  Ie = max(I);

          J = S.grid(dg).J_psi(ON); Js = min(J);  Je = max(J);

          Icorners = [Is Ie]; Jcorners = [Js Je];

          contact.corners.okay = true;
          contact.corners.Xdg  = G(dg).lon_psi(Icorners,Jcorners);
          contact.corners.Ydg  = G(dg).lat_psi(Icorners,Jcorners);
          contact.corners.Idg  = [Is Ie Ie Is];
          contact.corners.Jdg  = [Js Js Je Je];

% If any of edge boundary vectors are empty, the refined receiver grid
% is inside of more that one nested grid.  That is, the application has
% more than two nesting layers with telescoping refinement grids.
% For example, a refined grid (layer 3) is inside of another refined
% grid (layer 2), which in terms is inside of a coarser grid (layer 1).
% In this case the refined grid in layer 3 is not connected directly
% to the coarser grid in layer 1.  This will be an invalid contact
% region and we need to discard it.

          Jwest  = Js+1:Je-1;
          if (~isempty(Jwest))
            Iwest  = ones(size(Jwest )).*Is;
            B(iwest ).ind = sub2ind(size(S.grid(dg).I_psi),Iwest ,Jwest );
          else
            set_boundary = false;
          end

          Isouth = Is+1:Ie-1;
          if (~isempty(Isouth))
            Jsouth = ones(size(Isouth)).*Js;
            B(isouth).ind = sub2ind(size(S.grid(dg).I_psi),Isouth,Jsouth);
          else
            set_boundary = false;
          end

          Jeast  = Js+1:Je-1;
          if (~isempty(Jeast))
            Ieast  = ones(size(Jeast )).*Ie;
            B(ieast ).ind = sub2ind(size(S.grid(dg).I_psi),Ieast ,Jeast );
          else
            set_boundary = false;
          end

          Inorth = Is+1:Ie-1;
          if (~isempty(Inorth))
            Jnorth = ones(size(Inorth)).*Je;
            B(inorth).ind = sub2ind(size(S.grid(dg).I_psi),Inorth,Jnorth);
          else
            set_boundary = false;
          end

          if (set_boundary)
            for ib=1:4
              contact.boundary(ib).okay  = true;
              contact.boundary(ib).match = true(size(B(ib).ind));
              contact.boundary(ib).Xdg   = G(dg).lon_psi(B(ib).ind);
              contact.boundary(ib).Ydg   = G(dg).lat_psi(B(ib).ind);
              contact.boundary(ib).Idg   = S.grid(dg).I_psi(B(ib).ind);
              contact.boundary(ib).Jdg   = S.grid(dg).J_psi(B(ib).ind);
            end
            connected(dg,rg) = true;
            connected(rg,dg) = true;
          else
            for ib=1:4
              contact.boundary(ib).okay  = false;
              contact.boundary(ib).match = [];
              contact.boundary(ib).Xdg   = [];
              contact.boundary(ib).Ydg   = [];
              contact.boundary(ib).Idg   = [];
              contact.boundary(ib).Jdg   = [];
            end 
          end

        else

          [~,ON] = inpolygon(G(dg).x_psi, G(dg).y_psi,                  ...
                             S.grid(rg).perimeter.X_psi,                ...
                             S.grid(rg).perimeter.Y_psi);

          I = S.grid(dg).I_psi(ON); Is = min(I);  Ie = max(I);

          J = S.grid(dg).J_psi(ON); Js = min(J);  Je = max(J);

          Icorners = [Is Ie]; Jcorners = [Js Je];

          contact.corners.okay = true;
          contact.corners.Xdg  = G(dg).x_psi(Icorners,Jcorners);
          contact.corners.Ydg  = G(dg).y_psi(Icorners,Jcorners);
          contact.corners.Idg  = [Is Ie Ie Is];
          contact.corners.Jdg  = [Js Js Je Je];

% If any of edge boundary vectors are empty, the refined receiver grid
% is inside of more that one nested grid.  That is, the application has
% more than two nesting layers with telescoping refinement grids.
% For example, a refined grid (layer 3) is inside of another refined
% grid (layer 2), which in terms is inside of a coarser grid (layer 1).
% In this case the refined grid in layer 3 is not connected directly
% to the coarser grid in layer 1.  This will be an invalid contact
% region and we need to discard it.

          Jwest  = Js+1:Je-1;
          if (~isempty(Jwest))
            Iwest  = ones(size(Jwest )).*Is;
            B(iwest ).ind = sub2ind(size(S.grid(dg).I_psi),Iwest ,Jwest );
          else
            set_boundary = false;
          end

          Isouth = Is+1:Ie-1;
          if (~isempty(Isouth))
            Jsouth = ones(size(Isouth)).*Js;
            B(isouth).ind = sub2ind(size(S.grid(dg).I_psi),Isouth,Jsouth);
          else
            set_boundary = false;
          end

          Jeast  = Js+1:Je-1;
          if (~isempty(Jeast))
            Ieast  = ones(size(Jeast )).*Ie;
            B(ieast ).ind = sub2ind(size(S.grid(dg).I_psi),Ieast ,Jeast );
          else
            set_boundary = false;
          end

          Inorth = Is+1:Ie-1;
          if (~isempty(Inorth))
            Jnorth = ones(size(Inorth)).*Je;
            B(inorth).ind = sub2ind(size(S.grid(dg).I_psi),Inorth,Jnorth);
          else
            set_boundary = false;
          end

          if (set_boundary)
            for ib=1:4
              contact.boundary(ib).okay  = true;
              contact.boundary(ib).match = true(size(B(ib).ind));
              contact.boundary(ib).Xdg   = G(dg).x_psi(B(ib).ind);
              contact.boundary(ib).Ydg   = G(dg).y_psi(B(ib).ind);
              contact.boundary(ib).Idg   = S.grid(dg).I_psi(B(ib).ind);
              contact.boundary(ib).Jdg   = S.grid(dg).J_psi(B(ib).ind);
            end
            connected(dg,rg) = true;
            connected(rg,dg) = true;
          else
            for ib=1:4
              contact.boundary(ib).okay  = false;
              contact.boundary(ib).match = [];
              contact.boundary(ib).Xdg   = [];
              contact.boundary(ib).Ydg   = [];
              contact.boundary(ib).Idg   = [];
              contact.boundary(ib).Jdg   = [];
            end
          end
        
        end
      end

% Otherwise, determine if any of the donor grid boundary edges lay on the
% receiver grid perimeter.

      if (S.grid(rg).refine_factor == 0)
        for ib=1:4
          [~,ON] = inpolygon(S.grid(dg).boundary(ib).X,                 ...
                             S.grid(dg).boundary(ib).Y,                 ...
                             S.grid(rg).perimeter.X_psi,                ...
                             S.grid(rg).perimeter.Y_psi);
          if (any(ON))
            myindex = S.grid(dg).boundary(ib).index;
            contact.boundary(ib).okay  = true;
            contact.boundary(ib).match = [];
            contact.boundary(ib).index = myindex(ON);
            if (spherical)
              contact.boundary(ib).Xdg = G(dg).lon_psi(myindex(ON));
              contact.boundary(ib).Ydg = G(dg).lat_psi(myindex(ON));
            else
              contact.boundary(ib).Xdg = G(dg).x_psi(myindex(ON));
              contact.boundary(ib).Ydg = G(dg).y_psi(myindex(ON));
            end
            contact.boundary(ib).Idg = S.grid(dg).I_psi(myindex(ON));
            contact.boundary(ib).Jdg = S.grid(dg).J_psi(myindex(ON));
          else
            contact.boundary(ib).okay  = false;
            contact.boundary(ib).match = [];
            contact.boundary(ib).index = [];
            contact.boundary(ib).Xdg   = [];
            contact.boundary(ib).Ydg   = [];
            contact.boundary(ib).Idg   = [];
            contact.boundary(ib).Jdg   = [];
          end
        end
      end

% Determine if donor and receiver grids are coindident.  The logic below
% maybe inefficient but it is difficult to vectorize with uneven vectors.
% It matches the contact points at the boundary between donor and
% receiver grids.

      if (contact.corners.okay && ~contact.refinement)
        connected(dg,rg) = true;
        connected(rg,dg) = true;

        for ib=1:4
          ir = adjacent(ib);

          if (contact.boundary(ib).okay)
            dlength = length(S.grid(dg).boundary(ib).X);
            rlength = length(S.grid(rg).boundary(ir).X);

            dmatch  = false([1 dlength]);
            rmatch  = false([1 rlength]);

            icount = 0;     
            for n=1:dlength
              Xdg = S.grid(dg).boundary(ib).X(n);
              Ydg = S.grid(dg).boundary(ib).Y(n);
              for m=1:rlength
                Xrg = S.grid(rg).boundary(ir).X(m);
                Yrg = S.grid(rg).boundary(ir).Y(m);
                if ((Xdg == Xrg) && (Ydg == Yrg))
                  icount = icount+1;
                  dmatch(n) = true;
                  rmatch(m) = true;
                end
              end
            end

% There are coincident points.

            if (icount > 0)
              contact.coincident = true;
            end

% Either all the point along donor and receiver grids boundary edge
% are coincident (mosaic grids connectivity) or only few points are
% coincident (composite grids connectivity).

            if (dlength == icount && rlength == icount)
              contact.composite = true;
              contact.mosaic = true;
            else
              contact.composite = true;
            end

% Store indices of points in the boundary that are coincident.

            contact.boundary(ib).match = dmatch;
          else
            contact.boundary(ib).match = [];
          end
        end
      end

% Determine if grids are connected.

      if (contact.interior.okay)
        if (connected(dg,rg) || connected(rg,dg))
          load_connectivity = true;
        else
          load_connectivity = false;
        end
      else
        load_connectivity = false;
      end

% Load connectivity information into structure S if current grids
% 'dg' and 'rg' are connected.

      if (load_connectivity)
        cr = cr + 1;                          % contact region number
        S.contact(cr) = contact;
        disp([blanks(5), num2str(cr, '%2.2i'),                          ...
              blanks(8), num2str(dg, '%2.2i'),                          ...
              blanks(8), num2str(rg, '%2.2i')]);
      end
      clear contact

    end  % end of condional dg ~= rg
  end    % end of rg receiver grid loop
end      % end of dg donor grid loop

% Determine hybrid nesting application: coordinates should be
% identical in the domain in common. The substraction of the
% X-coordinate should be zero or less than epsilon.

if (S.Ngrids == 2)
  cr = 1;
  dg = S.contact(cr).donor_grid;
  rg = S.contact(cr).receiver_grid;
  if (G(dg).refine_factor <= 1 && G(rg).refine_factor <= 1)
    Imin(rg) = min(S.contact(cr).corners.Idg);
    Imax(rg) = max(S.contact(cr).corners.Idg);
    Jmin(rg) = min(S.contact(cr).corners.Jdg);
    Jmax(rg) = max(S.contact(cr).corners.Jdg);
    if (S.spherical)
      X=G(dg).lon_rho(Imin(rg):Imax(rg)+1, Jmin(rg):Jmax(rg)+1);
      Xmin=min(X(:)-G(rg).lon_rho(:));
      Xmax=max(X(:)-G(rg).lon_rho(:));
    else
      X=G(dg).x_rho(Imin(rg):Imax(rg)+1, Jmin(rg):Jmax(rg)+1);
%     Xmin=min(X(:)-G(rg).x_rho(:));  jcw orig
%     Xmax=max(X(:)-G(rg).x_rho(:));  jcw orig
      Xmin=min(min(X(:,:)-G(rg).x_rho(1:Imax(rg)+1-Imin(rg),:)));
      Xmax=max(max(X(:,:)-G(rg).x_rho(1:Imax(rg)+1-Imin(rg),:)));
    end
    if (Xmin < eps && Xmax < eps)
      S.contact(cr).hybrid   = true;
      S.contact(cr+1).hybrid = true;
    end
  end
end

disp(blanks(2));

return
