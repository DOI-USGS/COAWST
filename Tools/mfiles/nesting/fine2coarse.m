function C = fine2coarse(Ginp,Gout,Gfactor,varargin)

%
% FINE2COARSE:  Creates a coarser resolution ROMS Grid NetCDF file
%
% C = fine2coarse(Ginp,Gout,Gfactor,Imin,Imax,Jmin,Jmax)
%
% Given a fine resolution Grid NetCDF file (Ginp), this function creates
% a coarser resolution grid in the region specified by the finer grid
% coordinates (Imin,Jmin) and (Imax,Jmax). Notice that (Imin,Jmin), and
% (Imax,Jmax) indices are in terms of the PSI-points because it actually
% defines the physical boundaries of the coarser grid. The grid coarseness
% coefficient is specified with Gfactor.
%
% If Imin, Imax, Jmin, and Jmax is not provided, this function extracts
% the largest coarse grid possible within the finer grid.  If these
% PSI-indices are provided make sure that you pick the correct value
% from the following set for consistency between coarse and fine grids:
%
%    Imin, Imax      Any value from    Gfactor : Gfactor :L-Gfactor
%    Jmin, Jmax      Any value from    Gfactor : Gfactor :M-Gfactor
%
% where L and M are the number of PSI-points in the I- and J-directions.
%
% On Input:
%
%    Ginp       Input  finer   Grid NetCDF file name (string)
%    Gout       Output coarser Grid NetCDF file name (string)
%    Gfactor    Grid coarseness factor (3,5,7,9,11,13,15,...)
%    Imin       Finer grid lower-left  I-coordinate (PSI-point)
%    Imax       Finer grid upper-right I-coordinate (PSI-point)
%    Jmin       Finer grid lower-left  J-coordinate (PSI-point)
%    Jmax       Finer grid upper-right J-coordinate (PSI-point)
%
% On Output:
%
%    C          Coaser resolution Grid structure
%

% svn $Id: fine2coarse.m 729 2014-04-03 16:59:08Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
  
% If applicable, get larger grid structure.

if (~isstruct(Ginp)),
  F = get_roms_grid(Ginp);
  InpName = Ginp;
else
  F = Ginp;
  InpName = Ginp.grid_name;
end

% Set data sampling indices (PSI-points).

[Lp,Mp] = size(F.h);

L = Lp-1;
M = Mp-1;

% Set offset "half" factor which is used to extract boundary points
% outside of the PSI-points perimeter defined by Imin, Imax, Jmin,
% and Jmax.

half = floor(Gfactor-1)/2;

Imin = 1;                           % Default PSI-points range to apply
Imax = L-floor(mod(L,Gfactor)/2);   % coarseness to the full input grid
Jmin = 1;
Jmax = M-floor(mod(M,Gfactor)/2);

switch numel(varargin)
  case 1
    Imin = varargin{1};
  case 2
    Imin = varargin{1};
    Imax = varargin{2};
  case 3
    Imin = varargin{1};
    Imax = varargin{2};
    Jmin = varargin{3};
  case 4
    Imin = varargin{1};
    Imax = varargin{2};
    Jmin = varargin{3};
    Jmax = varargin{4};
end

% Check coarseness factor.

legal = abs([3:2:27] - Gfactor) < 4*eps;
if (~any(legal)),
  error([' FINE2COARSE: illegal coarseness factor, Gfactor = ',         ...
         num2str(Gfactor)]);
end

% Check extraction region.

if (Imin >= Imax),
  error([' FINE2COARSE: Imin >= Imax,    ',                             ...
         ' Imin = ', num2str(Imin),                                     ...
         ' Imax = ', num2str(Imax)]);
end

if (Jmin >= Jmax),
  error([' FINE2COARSE: Jmin >= Jmax,    ',                             ...
         ' Jmin = ', num2str(Jmin),                                     ...
         ' Jmax = ', num2str(Jmax)]);
end

if (Imin < 1),
  error([' FINE2COARSE: Imin < 1,   ',                                  ...
         ' Imin = ', num2str(Imax)]);
end

if (Imax > L),
  error([' FINE2COARSE: Imax > L,    ',                                 ...
         ' Imax = ', num2str(Imax),                                     ...
         ' L = ', num2str(L)]);
end

if (Jmin < 1),
  error([' FINE2COARSE: Jmin < 1,   ',                                  ...
         ' Jmin = ', num2str(Imax)]);
end

if (Jmax > M),
  error([' FINE2COARSE: Jmax > M,    ',                                 ...
         ' Jmax = ', num2str(Jmax),                                     ...
         ' M = ', num2str(M)]);
end

% Set grid variables to process.

grd_vars = {'h', 'f', 'angle', 'pm', 'pn'};

if (F.curvilinear),
  field_list = [grd_vars, 'dmde', 'dndx'];
end

grd_vars = [grd_vars, 'x_rho', 'y_rho', 'x_psi', 'y_psi',               ...
                      'x_u', 'y_u', 'x_v', 'y_v'};

if (F.spherical),
  grd_vars = [grd_vars, 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',     ...
                        'lon_u', 'lat_u', 'lon_v', 'lat_v'];
end

grd_vars = [grd_vars, 'mask_rho', 'mask_psi', 'mask_u', 'mask_v'];

grd_vars = [grd_vars, 'hraw'];               % needs to be last

%--------------------------------------------------------------------------
% Extract coaser grid.  Only valid for C-grid...
%--------------------------------------------------------------------------

% Set extraction ranges. The "half" offset value is to extract RHO-points
% boundaries (all four edges), U-points southern and northern boundary
% edges, and V-points western and eastern boundary edges. All these points
% are outside the fine grid perimeter defined by Imin, Imax, Jmin, and
% Jmax.

IminP = Imin;            ImaxP = Imax;
IminR = Imin-(half+1);   ImaxR = Imax+(half+1);
IminU = Imin;            ImaxU = Imax;
IminV = Imin-(half+1);   ImaxV = Imax+(half+1);

JminP = Jmin;            JmaxP = Jmax;
JminR = Jmin-(half+1);   JmaxR = Jmax+(half+1);
JminU = Jmin-(half+1);   JmaxU = Jmax+(half+1);
JminV = Jmin;            JmaxV = Jmax;

% Set fine and coarse grid ROMS I-indices to extract.

IindexP = false([1 L  ]);
IindexR = false([1 L+1]);
IindexU = false([1 L  ]);
IindexV = false([1 L+1]);

if (Imin == 1),
  IpC = half:Gfactor:L;
  IrC = 1   :Gfactor:L;
  IuC = half:Gfactor:L;
  IvC = 1   :Gfactor:L;
else
  IpC = Imin     :Gfactor:Imax;
  IrC = Imin-half:Gfactor:Imax;
  IuC = Imin     :Gfactor:Imax;
  IvC = Imin-half:Gfactor:Imax;
end

IindexP(IpC) = true;
IindexR(IrC) = true;
IindexU(IuC) = true;
IindexV(IvC) = true;

% Set fine and coarse grid ROMS I-indices to extract.

JindexP = false([1 M  ]);
JindexR = false([1 M+1]);
JindexU = false([1 M+1]);
JindexV = false([1 M  ]);

if (Jmin == 1),
  JpC = half:Gfactor:M;
  JrC = 1   :Gfactor:M;
  JuC = 1   :Gfactor:M;
  JvC = half:Gfactor:M;
else
  JpC = Jmin     :Gfactor:Jmax;
  JrC = Jmin-half:Gfactor:Jmax;
  JuC = Jmin-half:Gfactor:Jmax;
  JvC = Jmin     :Gfactor:Jmax;
end

JindexP(JpC) = true;
JindexR(JrC) = true;
JindexU(JuC) = true;
JindexV(JvC) = true;

% Initialize several subsomain structure paameters.

C.grid_name = Gout;
C.Lm = length(IpC)-1;
C.Mm = length(JpC)-1;
C.spherical = F.spherical;

if isfield(F, 'uniform'),
  C.uniform = F.uniform;
else
  C.uniform = 0;
end

% Extract grid variables.

for value = grd_vars,
  field = char(value);
  got.(field) = false;
  if (isfield(F,field)),
    switch (field)
      case {'lon_psi', 'lat_psi', 'mask_psi', 'x_psi', 'y_psi'}
        C.(field) = F.(field)(IindexP,JindexP);
        got.(field) = true;     
      case {'lon_u', 'lat_u', 'mask_u', 'x_u', 'y_u'}
        C.(field) = F.(field)(IindexU,JindexU);
        got.(field) = true;     
      case {'lon_v', 'lat_v', 'mask_v', 'x_v', 'y_v'}
        C.(field) = F.(field)(IindexV,JindexV);
          got.(field) = true;     
      case {'hraw'}     
        C.(field) = F.(field)(IindexR,JindexR,:);
        got.(field) = true;     
      otherwise
        C.(field) = F.(field)(IindexR,JindexR);
        got.(field) = true;     
    end
  end  
end

% Get grid lengths.

if (got.x_psi && got.y_psi),
  C.xl = max(C.x_psi(:)) - min(C.x_psi(:));
  C.el = max(C.y_psi(:)) - min(C.y_psi(:));
else
  C.xl = 0;
  C.el = 0;
end

% Recompute metrics at coarser resolution.  We cannot extract their
% values because grid spacing is larger by factor of Gfactors.
% The computation of metrics from discrete point is subject to
% roundoff.  There is not much that we can do here.  The roundoff
% is small and of the order 1.0E-16 (eps value).

disp(' ');
if (G.spherical),
  GreatCircle = true;
  disp('Computing grid spacing: great circle distances');
else
  GreatCircle = false;
  disp('Computing grid spacing: Cartesian distances');
end

[C.pm, C.pn, C.dndx, C.dmde]=grid_metrics(C, GreatCircle);

%--------------------------------------------------------------------------
% Create subdomain Grid NetCDF file and write out data.
%--------------------------------------------------------------------------

% Set number of grid points.

[Lp,Mp] = size(C.h);                 % RHO-points

disp(' ');
disp(['Number of points:',                                              ...
      ' Fine   Grid = ', num2str(F.Lm+2) ,' x ', num2str(F.Mm+2),       ...
      ',  Coarse Grid = ', num2str(C.Lm+2),' x ', num2str(C.Mm+2)]);

% Create ROMS Grid NetCDF file.

NewFile = true;

status = c_grid(Lp, Mp, Gout, NewFile, C.spherical);
if (status ~= 0), return, end

% Set global attributes.

status = nc_attadd(Gout, 'parent_grid', InpName);
if (status ~= 0), return, end

status = nc_attadd(Gout, 'parent_Imin', int32(Imin));
if (status ~= 0), return, end
  
status = nc_attadd(Gout, 'parent_Imax', int32(Imax));
if (status ~= 0), return, end

status = nc_attadd(Gout, 'parent_Jmin', int32(Jmin));
if (status ~= 0), return, end

status = nc_attadd(Gout, 'parent_Jmax', int32(Jmax));
if (status ~= 0), return, end

status = nc_attadd(Gout, 'coarse_factor', int32(Gfactor));
if (status ~= 0), return, end

status = nc_attdel(Gout, 'history');
if (status ~= 0), return, end

history = ['GRID file created using Matlab script: fine2coarse, ',      ...
           date_stamp];
status = nc_attadd(Gout, 'history', history);
if (status ~= 0), return, end

% Write out coarse resolution grid variables.

disp(['Writing subdomain grid variables into: ', Gout]);
disp(' ');

status = nc_write (Gout, 'spherical', C.spherical);
if (status ~= 0), return, end

status = nc_write (Gout, 'xl', C.xl);
if (status ~= 0), return, end

status = nc_write (Gout, 'el', C.el);
if (status ~= 0), return, end

if (got.hraw),
  bath = size(C.hraw,3);
  for rec=1:bath,
    status = nc_write (Gout, 'hraw', C.hraw, rec);
    if (status ~= 0), return, end
  end
else
  status = nc_write (Gout, 'hraw', C.h, 1);
  if (status ~= 0), return, end
end

for value = 1:length(grd_vars)-1,
  field = char(grd_vars(value));
  if (got.(field)),
    status = nc_write (Gout, field, C.(field));
    if (status ~= 0), return, end
  end  
end,

%--------------------------------------------------------------------------
% Add coastline data if available.
%--------------------------------------------------------------------------

if (isfield(F, 'lon_coast') && isfield(F, 'lat_coast')),
  add_coastline (Gout, F.lon_coast, F.lat_coast);
end

%--------------------------------------------------------------------------
% Get full extracted grid structure.
%--------------------------------------------------------------------------

C = get_roms_grid(Gout);

return
