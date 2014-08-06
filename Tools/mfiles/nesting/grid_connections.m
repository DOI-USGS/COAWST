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
%    S.contact(cr).mosaic                  - Mosaic grid switch
%    S.contact(cr).refinement              - Refinement grid switch
%
%    S.contact(cr).interior.okey           - true/false logical
%    S.contact(cr).interior.Xdg(:)         - (X,Y) coordinates and (I,J)
%    S.contact(cr).interior.Ydg(:)           indices of donor grid points
%    S.contact(cr).interior.Idg(:)           inside the receiver grid
%    S.contact(cr).interior.Jdg(:)           perimeter, [] if false 
%
%    S.contact(cr).corners.okey            - true/false logical
%    S.contact(cr).corners.Xdg(:)          - (X,Y) coordinates and (I,J)
%    S.contact(cr).corners.Ydg(:)            indices of donor grid points
%    S.contact(cr).corners.Idg(:)            corners laying on receiver
%    S.contact(cr).corners.Idg(:)            grid perimeter, [] if false
%
%    S.contact(cr).boundary(ib).okey       - true/false logical
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

% svn $Id: grid_connections.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
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

adjacent = [ieast, inorth, iwest, isouth];    % Receiver grid boundary
spherical = S.spherical;                      % spherical grid switch

%--------------------------------------------------------------------------
% Loop over all contact regions.
%--------------------------------------------------------------------------

cr = 0;

for dg=1:S.Ngrids,
  for rg=1:S.Ngrids,
    if (dg ~= rg),

      contact.donor_grid    = dg;             % donor grid number
      contact.receiver_grid = rg;             % receiver grid number
      contact.coincident    = false;
      contact.composite     = false;
      contact.mosaic        = false;

      if (S.grid(dg).refine_factor > 0 || S.grid(rg).refine_factor > 0),
        contact.refinement = true;
      else
        contact.refinement = false;
      end  

% Determine which points from donor grid are inside the reciever grid.
% If the receiver is a refinment grid, consider only donor RHO-points
% outside the receiver grid. Otherwise, consider only points inside
% the receiver grid.

      if (spherical),
        [IN,~] = inpolygon(G(dg).lon_rho, G(dg).lat_rho,                ...
                           S.grid(rg).perimeter.X_psi,                  ...
                           S.grid(rg).perimeter.Y_psi);
        if (any(IN(:))),
          X = G(dg).lon_rho(:);
          Y = G(dg).lat_rho(:);
          I = S.grid(dg).I_rho(:);
          J = S.grid(dg).J_rho(:);

          contact.interior.okey = true;
          if (S.grid(rg).refine_factor > 0),
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
          contact.interior.okey  = false;
          contact.interior.Xdg   = [];
          contact.interior.Ydg   = [];
          contact.interior.Idg   = [];
          contact.interior.Jdg   = [];
        end
      else
        [IN,~] = inpolygon(G(dg).x_rho, G(dg).y_rho,                    ...
                           S.grid(rg).perimeter.X_psi,                  ...
                           S.grid(rg).perimeter.Y_psi);
        if (any(IN(:))),
          X = G(dg).x_rho(:);
          Y = G(dg).y_rho(:);
          I = S.grid(dg).I_rho(:);
          J = S.grid(dg).J_rho(:);

          contact.interior.okey = true;
          if (S.grid(rg).refine_factor > 0),
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
          contact.interior.okey  = false;
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
      if (any(ON)),
        contact.corners.okey = true;
        myindex = S.grid(dg).corners.index;
        if (spherical),
          contact.corners.Xdg = G(dg).lon_psi(myindex(ON));
          contact.corners.Ydg = G(dg).lat_psi(myindex(ON));
        else
          contact.corners.Xdg = G(dg).x_psi(myindex(ON));
          contact.corners.Ydg = G(dg).y_psi(myindex(ON));
        end
        contact.corners.Idg   = S.grid(dg).I_psi(myindex(ON));
        contact.corners.Jdg   = S.grid(dg).J_psi(myindex(ON));
      else
        contact.corners.okey  = false;
        contact.corners.Xdg   = [];
        contact.corners.Ydg   = [];
        contact.corners.Idg   = [];
        contact.corners.Jdg   = [];
      end

% If receiver is a refinement grid from donor, determine which donor
% grid points lay in the receiver grid perimeter.

      if (S.grid(rg).refine_factor > 0 && contact.interior.okey),
        if (spherical)
          [~,ON] = inpolygon(G(dg).lon_psi, G(dg).lat_psi,              ...
                             S.grid(rg).perimeter.X_psi,                ...
                             S.grid(rg).perimeter.Y_psi);

          I = S.grid(dg).I_psi(ON); Is = min(I);  Ie = max(I);

          J = S.grid(dg).J_psi(ON); Js = min(J);  Je = max(J);

          Icorners = [Is Ie]; Jcorners = [Js Je];

          contact.corners.okey = true;
          contact.corners.Xdg  = G(dg).lon_psi(Icorners,Jcorners);
          contact.corners.Ydg  = G(dg).lat_psi(Icorners,Jcorners);
          contact.corners.Idg  = [Is Ie Ie Is];
          contact.corners.Jdg  = [Js Js Je Je];

          Jwest  = Js+1:Je-1;    Iwest  = ones(size(Jwest )).*Is;
          B(iwest ).ind = sub2ind(size(S.grid(dg).I_psi), Iwest , Jwest );

          Isouth = Is+1:Ie-1;    Jsouth = ones(size(Isouth)).*Js;
          B(isouth).ind = sub2ind(size(S.grid(dg).I_psi), Isouth, Jsouth);

          Jeast  = Js+1:Je-1;    Ieast  = ones(size(Jeast )).*Ie;
          B(ieast ).ind = sub2ind(size(S.grid(dg).I_psi), Ieast , Jeast );

          Inorth = Is+1:Ie-1;    Jnorth = ones(size(Inorth)).*Je;
          B(inorth).ind = sub2ind(size(S.grid(dg).I_psi), Inorth, Jnorth);
          
          for ib=1:4,
            contact.boundary(ib).okey  = true;
            contact.boundary(ib).match = true(size(B(ib).ind));
            contact.boundary(ib).Xdg   = G(dg).lon_psi(B(ib).ind);
            contact.boundary(ib).Ydg   = G(dg).lat_psi(B(ib).ind);
            contact.boundary(ib).Idg   = S.grid(dg).I_psi(B(ib).ind);
            contact.boundary(ib).Jdg   = S.grid(dg).J_psi(B(ib).ind);
          end
        else
          [~,ON] = inpolygon(G(dg).x_psi, G(dg).y_psi,                  ...
                             S.grid(rg).perimeter.X_psi,                ...
                             S.grid(rg).perimeter.Y_psi);

          I = S.grid(dg).I_psi(ON); Is = min(I);  Ie = max(I);

          J = S.grid(dg).J_psi(ON); Js = min(J);  Je = max(J);

          Icorners = [Is Ie]; Jcorners = [Js Je];

          contact.corners.okey = true;
          contact.corners.Xdg  = G(dg).x_psi(Icorners,Jcorners);
          contact.corners.Ydg  = G(dg).y_psi(Icorners,Jcorners);
          contact.corners.Idg  = [Is Ie Ie Is];
          contact.corners.Jdg  = [Js Js Je Je];

          Jwest  = Js+1:Je-1;    Iwest  = ones(size(Jwest )).*Is;
          B(iwest ).ind = sub2ind(size(S.grid(dg).I_psi), Iwest , Jwest );

          Isouth = Is+1:Ie-1;    Jsouth = ones(size(Isouth)).*Js;
          B(isouth).ind = sub2ind(size(S.grid(dg).I_psi), Isouth, Jsouth);

          Jeast  = Js+1:Je-1;    Ieast  = ones(size(Jeast )).*Ie;
          B(ieast ).ind = sub2ind(size(S.grid(dg).I_psi), Ieast , Jeast );

          Inorth = Is+1:Ie-1;    Jnorth = ones(size(Inorth)).*Je;
          B(inorth).ind = sub2ind(size(S.grid(dg).I_psi), Inorth, Jnorth);

          for ib=1:4,
            contact.boundary(ib).okey  = true;
            contact.boundary(ib).match = true(size(B(ib).ind));
            contact.boundary(ib).Xdg   = G(dg).x_psi(B(ib).ind);
            contact.boundary(ib).Ydg   = G(dg).y_psi(B(ib).ind);
            contact.boundary(ib).Idg   = S.grid(dg).I_psi(B(ib).ind);
            contact.boundary(ib).Jdg   = S.grid(dg).J_psi(B(ib).ind);
          end
        end
      end

% Otherwise, determine if any of the donor grid boundary edges lay on the
% receiver grid perimeter.

      if (S.grid(rg).refine_factor == 0),
        for ib=1:4,
          [~,ON] = inpolygon(S.grid(dg).boundary(ib).X,                 ...
                             S.grid(dg).boundary(ib).Y,                 ...
                             S.grid(rg).perimeter.X_psi,                ...
                             S.grid(rg).perimeter.Y_psi);
          if (any(ON)),
            myindex = S.grid(dg).boundary(ib).index;
            contact.boundary(ib).okey  = true;
            contact.boundary(ib).match = [];
            contact.boundary(ib).index = myindex(ON);
            if (spherical),
              contact.boundary(ib).Xdg = G(dg).lon_psi(myindex(ON));
              contact.boundary(ib).Ydg = G(dg).lat_psi(myindex(ON));
            else
              contact.boundary(ib).Xdg = G(dg).x_psi(myindex(ON));
              contact.boundary(ib).Ydg = G(dg).y_psi(myindex(ON));
            end
            contact.boundary(ib).Idg = S.grid(dg).I_psi(myindex(ON));
            contact.boundary(ib).Jdg = S.grid(dg).J_psi(myindex(ON));
          else
            contact.boundary(ib).okey  = false;
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

      if (contact.corners.okey && ~contact.refinement),
        for ib=1:4,
          ir = adjacent(ib);

          if (contact.boundary(ib).okey),
            dlength = length(S.grid(dg).boundary(ib).X);
            rlength = length(S.grid(rg).boundary(ir).X);

            dmatch  = false([1 dlength]);
            rmatch  = false([1 rlength]);

            icount = 0;     
            for n=1:dlength,
              Xdg = S.grid(dg).boundary(ib).X(n);
              Ydg = S.grid(dg).boundary(ib).Y(n);
              for m=1:rlength,
                Xrg = S.grid(rg).boundary(ir).X(m);
                Yrg = S.grid(rg).boundary(ir).Y(m);
                if ((Xdg == Xrg) && (Ydg == Yrg)),
                  icount = icount+1;
                  dmatch(n) = true;
                  rmatch(m) = true;
                end
              end
            end

% There are coincident points.

            if (icount > 0),
              contact.coincident = true;
            end

% Either all the point along donor and receiver grids boundary edge
% are coincident (mosaic grids connectivity) or only few points are
% coincident (composite grids connectivity).

            if (dlength == icount && rlength == icount),
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

% Load connectivity information into structure S if current grids
% 'dg' and 'rg' are connected.

      if (contact.interior.okey),
        cr = cr + 1;                          % contact region number
        S.contact(cr) = contact;
      end
      clear contact

    end  % end of condional dg ~= rg
  end    % end of rg receiver grid loop
end      % end of dg donor grid loop

return
