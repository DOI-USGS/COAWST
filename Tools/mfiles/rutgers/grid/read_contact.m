function [S] = read_contact(ncname)

%
% READ_CONTACT:  Reads ROMS Nested Grids Contact Points NetCDF file
%
% S = read_contact(ncname)
%
% This function reads in the Nested Grids Contact Point NetCDF file
% and loads data into a structure.
%
% On Input:
%
%    ncname      Contact Points NetCDF file name (string)
%
% On Output:
%
%    S           Nested grids Contact Points structure (struct array)
%

% svn $Id: read_contact.m 738 2014-10-14 21:49:14Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Inquire about the contents of NetCDF.
  
I = nc_inq(ncname);

% Get dimension parameters. 

Ngrids    = I.Dimensions(strcmp({I.Dimensions.Name},'Ngrids'   )).Length;
Ncontact  = I.Dimensions(strcmp({I.Dimensions.Name},'Ncontact' )).Length;
nLweights = I.Dimensions(strcmp({I.Dimensions.Name},'nLweights')).Length;
nQweights = I.Dimensions(strcmp({I.Dimensions.Name},'nQweights')).Length;
Ndatum    = I.Dimensions(strcmp({I.Dimensions.Name},'datum'    )).Length;

%--------------------------------------------------------------------------
% Check if Nested Grid NetCDF are available and get their grid structures.
%--------------------------------------------------------------------------

% Get nested grid NetCDF file names.

filenames = I.Attributes(strcmp({I.Attributes.Name},'grid_files')).Value;

% Remove newline control characters, sprintf('\n').

files = filenames(isstrprop(filenames, 'graphic'));  

ind = strfind(files, '.nc')+3;

Gnames  = {};
gotfile = false([1 Ngrids]);

ie = 0;
for ng=1:Ngrids
  is = ie+1;
  ie = ind(ng)-1;
  Gnames = [Gnames files(is:ie)];
  gotfile(ng) = exist(char(Gnames(ng)), 'file');
end

% Get nested grids structures. If the grid structure have the parent
% fields, remove them to have an array of similar structures.

if (all(gotfile)),
  parent = {'parent_grid',                                              ...
            'parent_Imin', 'parent_Imax',                               ...
            'parent_Jmin', 'parent_Jmax'};
  for ng=1:Ngrids,
    g = get_roms_grid(char(Gnames(ng)));
    if (isfield(g, 'parent_grid')),
      G(ng) = rmfield(g, parent);
    else
      G(ng) = g;
    end
  end
end

% Initialize structure if perimeter, boundary edges, and connectivity
% information.

if (all(gotfile)),
  S = grid_perimeter(G);
  S = grid_connections(G, S);
end

%--------------------------------------------------------------------------
% Build Contact points NetCDF structure.
%--------------------------------------------------------------------------

% Read in several grid parameters.

par_list = {'Lm', 'Mm',                                                 ...
            'coincident', 'composite', 'mosaic', 'refinement',          ...
            'refine_factor',                                            ...
            'donor_grid', 'receiver_grid',                              ...
            'NstrR', 'NendR',                                           ...
            'NstrU', 'NendU',                                           ...
            'NstrV', 'NendV'};

for value = par_list,
  field = char(value);
  P.(field) = ncread(ncname, field);
end

P.coincident = logical(P.coincident);
P.composite  = logical(P.composite);
P.mosaic     = logical(P.mosaic);
P.refinement = logical(P.refinement);

% Load values into structure.

if (all(gotfile)),
  S.Ndatum = Ndatum;
else
  S.Ngrids   = Ngrids;               % Grid NetCDF filenames in global
  S.Ncontact = Ncontact;             % attributes are not available.
  S.Nweights = Nweights;             % Therefore, build the S.grid
  S.Ndatum   = Ndatum;               % substructure with the few
                                     % informations that it is
  S.western_edge  = 1;               % available in the NetCDF file.
  S.southern_edge = 2;
  S.eastern_edge  = 3;
  S.northern_edge = 4;

  S.spherical = ncread(ncname, 'spherical');

  for ng=1:Ngrids
    S.grid(ng).filename = char(Gnames(ng));

    S.grid(ng).Lp = P.Lm(ng)+2;
    S.grid(ng).Mp = P.Mm(ng)+2;
  
    S.grid(ng).L  = P.Lm(ng)+1;
    S.grid(ng).M  = P.Mm(ng)+1;
  
    S.grid(ng).refine_factor = P.refine_factor(ng);
  end
end

% Read in contact point variables.

var_list = {'contact_region',                                           ...
            'on_boundary',                                              ...
            'Idg', 'Jdg',                                               ...
            'Irg', 'Jrg',                                               ...
            'Lweight', 'Qweight',                                       ...
            'h', 'f', 'pm', 'pn', 'dndx', 'dmde',                       ...
            'Xrg', 'Yrg', 'angle', 'mask'};

for value = var_list,
  field = char(value);
  C.(field) = double(ncread(ncname, field));
end

C.on_boundary = logical(C.on_boundary);

% Load values into structure.

for n=1:Ncontact,
  S.contact(n).donor_grid         = P.donor_grid(n);
  S.contact(n).receiver_grid      = P.receiver_grid(n);
  S.contact(n).coincident         = P.coincident(n);
  S.contact(n).composite          = P.composite(n);
  S.contact(n).mosaic             = P.mosaic(n);
  S.contact(n).refinement         = P.refinement(n);
  
  S.contact(n).point.Xrg_rho      = C.Xrg(P.NstrR(n):P.NendR(n));
  S.contact(n).point.Yrg_rho      = C.Yrg(P.NstrR(n):P.NendR(n));
  S.contact(n).point.Irg_rho      = C.Irg(P.NstrR(n):P.NendR(n));
  S.contact(n).point.Jrg_rho      = C.Jrg(P.NstrR(n):P.NendR(n));
  S.contact(n).point.Idg_rho      = C.Idg(P.NstrR(n):P.NendR(n));
  S.contact(n).point.Jdg_rho      = C.Jdg(P.NstrR(n):P.NendR(n));

  S.contact(n).point.Xrg_u        = C.Xrg(P.NstrU(n):P.NendU(n));
  S.contact(n).point.Yrg_u        = C.Yrg(P.NstrU(n):P.NendU(n));
  S.contact(n).point.Irg_u        = C.Irg(P.NstrU(n):P.NendU(n));
  S.contact(n).point.Jrg_u        = C.Jrg(P.NstrU(n):P.NendU(n));
  S.contact(n).point.Idg_u        = C.Idg(P.NstrU(n):P.NendU(n));
  S.contact(n).point.Jdg_u        = C.Jdg(P.NstrU(n):P.NendU(n));

  S.contact(n).point.Xrg_v        = C.Xrg(P.NstrV(n):P.NendV(n));
  S.contact(n).point.Yrg_v        = C.Yrg(P.NstrV(n):P.NendV(n));
  S.contact(n).point.Irg_v        = C.Irg(P.NstrV(n):P.NendV(n));
  S.contact(n).point.Jrg_v        = C.Jrg(P.NstrV(n):P.NendV(n));
  S.contact(n).point.Idg_v        = C.Idg(P.NstrV(n):P.NendV(n));
  S.contact(n).point.Jdg_v        = C.Jdg(P.NstrV(n):P.NendV(n));

  S.contact(n).point.boundary_rho = C.on_boundary(P.NstrR(n):P.NendR(n));
  S.contact(n).point.boundary_u   = C.on_boundary(P.NstrU(n):P.NendU(n));
  S.contact(n).point.boundary_v   = C.on_boundary(P.NstrV(n):P.NendV(n));

  S.contact(n).point.angle        = C.angle(P.NstrR(n):P.NendR(n));
  S.contact(n).point.f            = C.f    (P.NstrR(n):P.NendR(n));
  S.contact(n).point.h            = C.h    (P.NstrR(n):P.NendR(n));
  S.contact(n).point.pm           = C.pm   (P.NstrR(n):P.NendR(n));
  S.contact(n).point.pn           = C.pn   (P.NstrR(n):P.NendR(n));
  S.contact(n).point.dndx         = C.dndx (P.NstrR(n):P.NendR(n));
  S.contact(n).point.dmde         = C.dmde (P.NstrR(n):P.NendR(n));
  S.contact(n).point.mask_rho     = C.mask (P.NstrR(n):P.NendR(n));
  S.contact(n).point.mask_u       = C.mask (P.NstrU(n):P.NendU(n));
  S.contact(n).point.mask_v       = C.mask (P.NstrV(n):P.NendV(n));

  S.Lweights(n).H_rho             = C.Lweight(1:nLweights,              ...
                                              P.NstrR(n):P.NendR(n));
  S.Lweights(n).H_u               = C.Lweight(1:nLweights,              ...
                                              P.NstrU(n):P.NendU(n));
  S.Lweights(n).H_v               = C.Lweight(1:nLweights,              ...
                                              P.NstrV(n):P.NendV(n));

  S.Qweights(n).H_rho             = C.Qweight(1:nQweights,              ...
                                              P.NstrR(n):P.NendR(n));
  S.Qweights(n).H_u               = C.Qweight(1:nQweights,              ...
                                              P.NstrU(n):P.NendU(n));
  S.Qweights(n).H_v               = C.Qweight(1:nQweights,              ...
                                              P.NstrV(n):P.NendV(n));
end

return
