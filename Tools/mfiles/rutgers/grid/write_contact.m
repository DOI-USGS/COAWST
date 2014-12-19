function write_contact(ncname, S, G, varargin)

%
% WRITE_CONTACT:  Writes ROMS Nested Grids Contact Points to a NetCDF file
%
% write_contact(ncname, S, Lcreate)
%
% This function writes out the Nested Grids Contact Point data generated
% by either "contact.m" or "read_contact.m" into specified NetCDF file.
%
% On Input:
%
%    ncname      Contact Point NetCDF file name (string)
%
%    S           Nested grids Contact Points structure (struct array)
%                  by "contact.m" or "read_contact.m" 
%
%    G           Information grids structure (1 x Ngrids struct array)
%
%    Lcreate     Switch to create a new NetCDF file (optional,
%                                                    default true)
%

% svn $Id: write_contact.m 746 2014-12-15 23:24:27Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

Lcreate = true;

switch numel(varargin)
  case 1
    Lcreate = varargin{1};
end

Ngrids    = S.Ngrids;           % number of nested grids
Ncontact  = S.Ncontact;         % number of Contact Regions
nLweights = S.nLweights;        % Number of linear    interpolation weights
nQweights = S.nQweights;        % Number of quadratic interpolation weights
Ndatum    = S.Ndatum;           % total number of Contact Points
spherical = S.spherical;

FillValue = 1.0d+37;            % ROMS FillValue

newline   = sprintf('\n');      % newline in global attributes

for ng = 1:Ngrids,
  AreaAvg(ng) = mean(mean((1./G(ng).pm) .* (1./G(ng).pn)));
end  

%--------------------------------------------------------------------------
% Create Contact Points NetCDF file.
%--------------------------------------------------------------------------

if (Lcreate),

% Create NetCDF file.
  
  c_contact(ncname, spherical, Ngrids, Ndatum); 

% Define few global attributes.

  Gnames = [];
  for n=1:Ngrids,
    Gnames = [Gnames  blanks(1) newline S.grid(n).filename];
  end
  ncwriteatt(ncname, '/', 'grid_files', Gnames);

% History attribute.

  history = ['Contact points created by contact.m ',                    ...
             'and written by write_contact.m, ', date_stamp]; 
  ncwriteatt(ncname, '/', 'history', history);

end

%--------------------------------------------------------------------------
% Write out Contact Points data.
%--------------------------------------------------------------------------

% Spherical switch.

ncwrite(ncname, 'spherical', int32(S.spherical));

% Number of interior RHO-points.

ncwrite(ncname, 'Lm', int32([S.grid.Lp]-2));
ncwrite(ncname, 'Mm', int32([S.grid.Mp]-2));

% Information variables.

ncwrite(ncname, 'coincident', int32([S.contact.coincident]));
ncwrite(ncname, 'composite', int32([S.contact.composite]));
ncwrite(ncname, 'mosaic', int32([S.contact.mosaic]));
ncwrite(ncname, 'refinement', int32([S.contact.refinement] > 0));
ncwrite(ncname, 'refine_factor', int32([S.grid.refine_factor]));

% Contact Points vertical interpolation switch.

ncwrite(ncname, 'interpolate', int32(ones([1 Ncontact])));

% Contact region donor and receiver grid.

ncwrite(ncname, 'donor_grid', int32([S.contact.donor_grid]));
ncwrite(ncname, 'receiver_grid', int32([S.contact.receiver_grid]));

% Refinement grid extraction coordinates, if any.

spv  = -999;
Imin = ones([1 Ngrids]) .* spv;
Jmin = ones([1 Ngrids]) .* spv;
Imax = ones([1 Ngrids]) .* spv;
Jmax = ones([1 Ngrids]) .* spv;

for cr=1:Ncontact
  dg = S.contact(cr).donor_grid;
  rg = S.contact(cr).receiver_grid;
  if ((S.grid(rg).refine_factor > 0) && AreaAvg(dg) > AreaAvg(rg)),
    Imin(rg) = min(S.contact(cr).corners.Idg);
    Jmin(rg) = min(S.contact(cr).corners.Jdg);
    Imax(rg) = max(S.contact(cr).corners.Idg);
    Jmax(rg) = max(S.contact(cr).corners.Jdg);
  end
end

ncwrite(ncname, 'I_left',   int32(Imin));
ncwrite(ncname, 'I_right',  int32(Imax));
ncwrite(ncname, 'J_bottom', int32(Jmin));
ncwrite(ncname, 'J_top',    int32(Jmax));

% Number of contact point is each contact region.

NpointsR = zeros([1 Ncontact]);
NpointsU = zeros([1 Ncontact]);
NpointsV = zeros([1 Ncontact]);

NstrR = zeros([1 Ncontact]);
NendR = zeros([1 Ncontact]);

NstrU = zeros([1 Ncontact]);
NendU = zeros([1 Ncontact]);

NstrV = zeros([1 Ncontact]);
NendV = zeros([1 Ncontact]);

ic = 0;

for cr=1:Ncontact
  NpointsR(cr) = length(S.contact(cr).point.Irg_rho);
  NpointsU(cr) = length(S.contact(cr).point.Irg_u);
  NpointsV(cr) = length(S.contact(cr).point.Irg_v);

  if (cr == 1),
    NstrR(cr) = 1;
  else
    NstrR(cr) = NendV(cr-1) + 1;
  end
  NendR(cr) = NstrR(cr) + NpointsR(cr) - 1;
  NstrU(cr) = NendR(cr) + 1;
  NendU(cr) = NstrU(cr) + NpointsU(cr) - 1;
  NstrV(cr) = NendU(cr) + 1;
  NendV(cr) = NstrV(cr) + NpointsV(cr) - 1;
end

ncwrite(ncname, 'NstrR', int32(NstrR));
ncwrite(ncname, 'NendR', int32(NendR));

ncwrite(ncname, 'NstrU', int32(NstrU));
ncwrite(ncname, 'NendU', int32(NendU));

ncwrite(ncname, 'NstrV', int32(NstrV));
ncwrite(ncname, 'NendV', int32(NendV));

% Contact region for each contact point.

contact_region = zeros([1 Ndatum]);

for cr=1:Ncontact,
  contact_region(NstrR(cr):NendV(cr)) = cr;
end

ncwrite(ncname, 'contact_region', int32(contact_region));

% Donor grid cell indices.

Idg = NaN([1 Ndatum]);
Jdg = NaN([1 Ndatum]);

for cr=1:Ncontact,
  Idg(NstrR(cr):NendR(cr)) = S.contact(cr).point.Idg_rho;
  Idg(NstrU(cr):NendU(cr)) = S.contact(cr).point.Idg_u;
  Idg(NstrV(cr):NendV(cr)) = S.contact(cr).point.Idg_v;

  Jdg(NstrR(cr):NendR(cr)) = S.contact(cr).point.Jdg_rho;
  Jdg(NstrU(cr):NendU(cr)) = S.contact(cr).point.Jdg_u;
  Jdg(NstrV(cr):NendV(cr)) = S.contact(cr).point.Jdg_v;
end

ncwrite(ncname, 'Idg', int32(Idg));
ncwrite(ncname, 'Jdg', int32(Jdg));

% Receiver grid indices.

Irg = NaN([1 Ndatum]);
Jrg = NaN([1 Ndatum]);

for cr=1:Ncontact,
  Irg(NstrR(cr):NendR(cr)) = S.contact(cr).point.Irg_rho;
  Irg(NstrU(cr):NendU(cr)) = S.contact(cr).point.Irg_u;
  Irg(NstrV(cr):NendV(cr)) = S.contact(cr).point.Irg_v;

  Jrg(NstrR(cr):NendR(cr)) = S.contact(cr).point.Jrg_rho;
  Jrg(NstrU(cr):NendU(cr)) = S.contact(cr).point.Jrg_u;
  Jrg(NstrV(cr):NendV(cr)) = S.contact(cr).point.Jrg_v;
end

ncwrite(ncname, 'Irg', int32(Irg));
ncwrite(ncname, 'Jrg', int32(Jrg));

% Donor grid cell bottom-left (XI,ETA) coordinates.

xi_dg  = NaN([1 Ndatum]);
eta_dg = NaN([1 Ndatum]);

for cr=1:Ncontact,
  dg = S.contact(cr).donor_grid;
  [Ir,Jr] = size(S.grid(dg).I_rho);
  [Iu,Ju] = size(S.grid(dg).I_u);
  [Iv,Jv] = size(S.grid(dg).I_v);
  Rindex  = sub2ind([Ir, Jr],                                           ...
                      S.contact(cr).point.Idg_rho+1,                    ...
                      S.contact(cr).point.Jdg_rho+1);
  Uindex  = sub2ind([Iu, Ju],                                           ...
                      S.contact(cr).point.Idg_u,                        ...
                      S.contact(cr).point.Jdg_u+1);
  Vindex  = sub2ind([Iv, Jv],                                           ...
                      S.contact(cr).point.Idg_v+1,                      ...
                      S.contact(cr).point.Jdg_v);

  xi_dg(NstrR(cr):NendR(cr)) = S.grid(dg).XI_rho(Rindex);
  xi_dg(NstrU(cr):NendU(cr)) = S.grid(dg).XI_u  (Uindex);
  xi_dg(NstrV(cr):NendV(cr)) = S.grid(dg).XI_v  (Vindex);
  
  eta_dg(NstrR(cr):NendR(cr)) = S.grid(dg).ETA_rho(Rindex);
  eta_dg(NstrU(cr):NendU(cr)) = S.grid(dg).ETA_u  (Uindex);
  eta_dg(NstrV(cr):NendV(cr)) = S.grid(dg).ETA_v  (Vindex);
end

ncwrite(ncname, 'xi_dg',  xi_dg);
ncwrite(ncname, 'eta_dg', eta_dg);

% Receiver grid contact point (XI,ETA) coordinates.

xi_rg  = NaN([1 Ndatum]);
eta_rg = NaN([1 Ndatum]);

for cr=1:Ncontact,
  xi_rg(NstrR(cr):NendR(cr)) = S.contact(cr).point.xrg_rho;
  xi_rg(NstrU(cr):NendU(cr)) = S.contact(cr).point.xrg_u;
  xi_rg(NstrV(cr):NendV(cr)) = S.contact(cr).point.xrg_v;

  eta_rg(NstrR(cr):NendR(cr)) = S.contact(cr).point.erg_rho;
  eta_rg(NstrU(cr):NendU(cr)) = S.contact(cr).point.erg_u;
  eta_rg(NstrV(cr):NendV(cr)) = S.contact(cr).point.erg_v;
end

ncwrite(ncname, 'xi_rg',  xi_rg);
ncwrite(ncname, 'eta_rg', eta_rg);

% Donor grid cell bottom-left coordinates.

Xdg = NaN([1 Ndatum]);
Ydg = NaN([1 Ndatum]);

for cr=1:Ncontact,
  dg = S.contact(cr).donor_grid;
  [Ir,Jr] = size(S.grid(dg).I_rho);
  [Iu,Ju] = size(S.grid(dg).I_u);
  [Iv,Jv] = size(S.grid(dg).I_v);
  Rindex  = sub2ind([Ir, Jr],                                           ...
                      S.contact(cr).point.Idg_rho+1,                    ...
                      S.contact(cr).point.Jdg_rho+1);
  Uindex  = sub2ind([Iu, Ju],                                           ...
                      S.contact(cr).point.Idg_u,                        ...
                      S.contact(cr).point.Jdg_u+1);
  Vindex  = sub2ind([Iv, Jv],                                           ...
                      S.contact(cr).point.Idg_v+1,                      ...
                      S.contact(cr).point.Jdg_v);
  if (S.spherical),
    Xdg(NstrR(cr):NendR(cr)) = G(dg).lon_rho(Rindex);
    Xdg(NstrU(cr):NendU(cr)) = G(dg).lon_u  (Uindex);
    Xdg(NstrV(cr):NendV(cr)) = G(dg).lon_v  (Vindex);
  
    Ydg(NstrR(cr):NendR(cr)) = G(dg).lat_rho(Rindex);
    Ydg(NstrU(cr):NendU(cr)) = G(dg).lat_u  (Uindex);
    Ydg(NstrV(cr):NendV(cr)) = G(dg).lat_v  (Vindex);
  else
    Xdg(NstrR(cr):NendR(cr)) = G(dg).x_rho(Rindex);
    Xdg(NstrU(cr):NendU(cr)) = G(dg).x_u  (Uindex);
    Xdg(NstrV(cr):NendV(cr)) = G(dg).x_v  (Vindex);
  
    Ydg(NstrR(cr):NendR(cr)) = G(dg).y_rho(Rindex);
    Ydg(NstrU(cr):NendU(cr)) = G(dg).y_u  (Uindex);
    Ydg(NstrV(cr):NendV(cr)) = G(dg).y_v  (Vindex);
  end  
end

ncwrite(ncname, 'Xdg', Xdg);
ncwrite(ncname, 'Ydg', Ydg);

% Receiver grid contact points locations.

Xrg = NaN([1 Ndatum]);
Yrg = NaN([1 Ndatum]);

for cr=1:Ncontact,
  Xrg(NstrR(cr):NendR(cr)) = S.contact(cr).point.Xrg_rho;
  Xrg(NstrU(cr):NendU(cr)) = S.contact(cr).point.Xrg_u;
  Xrg(NstrV(cr):NendV(cr)) = S.contact(cr).point.Xrg_v;

  Yrg(NstrR(cr):NendR(cr)) = S.contact(cr).point.Yrg_rho;
  Yrg(NstrU(cr):NendU(cr)) = S.contact(cr).point.Yrg_u;
  Yrg(NstrV(cr):NendV(cr)) = S.contact(cr).point.Yrg_v;
end

ncwrite(ncname, 'Xrg', Xrg);
ncwrite(ncname, 'Yrg', Yrg);

% Horizontal linear interpolation weights.

Lweight = NaN([nLweights Ndatum]);

for cr=1:Ncontact,
  Lweight(1:nLweights, NstrR(cr):NendR(cr)) = S.Lweights(cr).H_rho;
  Lweight(1:nLweights, NstrU(cr):NendU(cr)) = S.Lweights(cr).H_u;
  Lweight(1:nLweights, NstrV(cr):NendV(cr)) = S.Lweights(cr).H_v;
end

ncwrite(ncname, 'Lweight', Lweight, [1 1]);

% Horizontal quadratic interpolation weights.

Qweight = NaN([nQweights Ndatum]);

for cr=1:Ncontact,
  Qweight(1:nQweights, NstrR(cr):NendR(cr)) = S.Qweights(cr).H_rho;
  Qweight(1:nQweights, NstrU(cr):NendU(cr)) = S.Qweights(cr).H_u;
  Qweight(1:nQweights, NstrV(cr):NendV(cr)) = S.Qweights(cr).H_v;
end

ncwrite(ncname, 'Qweight', Qweight, [1 1]);

% Several grid variables.

angle(1:Ndatum) = FillValue;
f    (1:Ndatum) = FillValue;
h    (1:Ndatum) = FillValue;
pm   (1:Ndatum) = FillValue;
pn   (1:Ndatum) = FillValue;
dndx (1:Ndatum) = FillValue;
dmde (1:Ndatum) = FillValue;

for cr=1:Ncontact,
  angle(NstrR(cr):NendR(cr)) = S.contact(cr).point.angle;
  f    (NstrR(cr):NendR(cr)) = S.contact(cr).point.f;
  h    (NstrR(cr):NendR(cr)) = S.contact(cr).point.h;
  pm   (NstrR(cr):NendR(cr)) = S.contact(cr).point.pm;
  pn   (NstrR(cr):NendR(cr)) = S.contact(cr).point.pn;
  dndx (NstrR(cr):NendR(cr)) = S.contact(cr).point.dndx;
  dmde (NstrR(cr):NendR(cr)) = S.contact(cr).point.dmde;
end

ncwrite(ncname, 'angle', angle);
ncwrite(ncname, 'f'    , f    );
ncwrite(ncname, 'h'    , h    );
ncwrite(ncname, 'pm'   , pm   );
ncwrite(ncname, 'pn'   , pn   );
ncwrite(ncname, 'dndx' , dndx );
ncwrite(ncname, 'dmde' , dmde );

% Land/Sea masking.

mask = NaN([1 Ndatum]);

for cr=1:Ncontact,
  mask(NstrR(cr):NendR(cr)) = S.contact(cr).point.mask_rho;
  mask(NstrU(cr):NendU(cr)) = S.contact(cr).point.mask_u;
  mask(NstrV(cr):NendV(cr)) = S.contact(cr).point.mask_v;
end

ncwrite(ncname, 'mask', mask);

%--------------------------------------------------------------------------
%  Determine which contact points are on receiver grid boundaries.
%--------------------------------------------------------------------------
%
% Contact points on receiver grid boundaries. The physical boundary
% is either at U-points (west and east edges) or V-points (south
% and north edges).

boundary = false([1 Ndatum]);

for cr=1:Ncontact,
  boundary(NstrU(cr):NendU(cr)) = S.contact(cr).point.boundary_u;
  boundary(NstrV(cr):NendV(cr)) = S.contact(cr).point.boundary_v;
end

on_boundary = zeros([1 Ndatum]);

for i=1:Ndatum,
  if (boundary(i)),
    cr=contact_region(i);
    rg=S.contact(cr).receiver_grid;
    for ib=1:4,
      xb = S.grid(rg).boundary(ib).Xuv;
      yb = S.grid(rg).boundary(ib).Yuv;
      ind = (abs(xb-Xrg(i)) < 4*eps & abs(yb-Yrg(i)) < 4*eps);
      if (any(ind)),
        on_boundary(i) = ib;
        break;
      end
    end
  end
end

ncwrite(ncname, 'on_boundary', int32([on_boundary]));

return
