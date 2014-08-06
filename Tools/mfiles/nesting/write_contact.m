function write_contact(ncname, S, varargin)

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
%    Lcreate     Switch to create a new NetCDF file (optional,
%                                                    default true)
%

% svn $Id: write_contact.m 711 2014-01-23 20:36:13Z arango $
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
Nweights  = S.Nweights;         % Number of horizontal weights
Ndatum    = S.Ndatum;           % total number of Contact Points
spherical = S.spherical;

FillValue = 1.0d+37;            % ROMS FillValue

newline   = sprintf('\n');      % newline in global attributes

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

% Horizontal interpolation weights.

Hweight = NaN([Nweights Ndatum]);

for cr=1:Ncontact,
  Hweight(1:Nweights, NstrR(cr):NendR(cr)) = S.weights(cr).H_rho;
  Hweight(1:Nweights, NstrU(cr):NendU(cr)) = S.weights(cr).H_u;
  Hweight(1:Nweights, NstrV(cr):NendV(cr)) = S.weights(cr).H_v;
end

ncwrite(ncname, 'Hweight', Hweight, [1 1]);

% Contact point locations.

X = NaN([1 Ndatum]);
Y = NaN([1 Ndatum]);

for cr=1:Ncontact,
  X(NstrR(cr):NendR(cr)) = S.contact(cr).point.Xrg_rho;
  X(NstrU(cr):NendU(cr)) = S.contact(cr).point.Xrg_u;
  X(NstrV(cr):NendV(cr)) = S.contact(cr).point.Xrg_v;

  Y(NstrR(cr):NendR(cr)) = S.contact(cr).point.Yrg_rho;
  Y(NstrU(cr):NendU(cr)) = S.contact(cr).point.Yrg_u;
  Y(NstrV(cr):NendV(cr)) = S.contact(cr).point.Yrg_v;
end

ncwrite(ncname, 'Xrg', X);
ncwrite(ncname, 'Yrg', Y);

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
      ind = (abs(xb-X(i)) < 4*eps & abs(yb-Y(i)) < 4*eps);
      if (any(ind)),
        on_boundary(i) = ib;
        break;
      end
    end
  end
end

ncwrite(ncname, 'on_boundary', int32([on_boundary]));

return
