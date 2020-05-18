function C = fine2coarse(Ginp, Gout, Gfactor, Lplot, varargin)

%
% FINE2COARSE:  Creates a coarser resolution ROMS Grid NetCDF file
%
% C = fine2coarse(Ginp, Gout, Gfactor, Lplot, Imin, Imax, Jmin, Jmax)
%
% Given a fine resolution Grid NetCDF file (Ginp), this function creates
% a coarser resolution grid in the region specified by the finer grid
% coordinates (Imin,Jmin) and (Imax,Jmax). Notice that (Imin,Jmin), and
% (Imax,Jmax) indices are in terms of the PSI-points because it actually
% defines the physical boundaries of the coarser grid. The grid coarseness
% coefficient is specified with Gfactor.
%
% If Imin, Imax, Jmin, and Jmax are not provided, this function extracts
% the largest coarse grid possible from the finer PSI-grid:
%
%    Imin = 1+half;
%    Imax = L-half+1;
%    Jmin = 1+half;
%    Jmax = M-half+1;
%
% where
%
%    half = (Gfactor-1)/2;
%
% here L and M are the finer grid number of PSI-points in the I- and
% J-directions, respectible.
%
%
% On Input:
%
%    Ginp       Input  finer   Grid NetCDF file name (string)
%
%    Gout       Output coarser Grid NetCDF file name (string)
%
%    Gfactor    Grid coarseness factor (3,5,7,9,11,13,15,...)
%
%    Lplot      Draw diagrams of fine and coarse grids.
%
%    Imin       Finer grid lower-left  I-coordinate (PSI-point)
%
%    Imax       Finer grid upper-right I-coordinate (PSI-point)
%
%    Jmin       Finer grid lower-left  J-coordinate (PSI-point)
%
%    Jmax       Finer grid upper-right J-coordinate (PSI-point)
%
% On Output:
%
%    C          Coaser resolution Grid structure
%

% svn $Id: fine2coarse.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
  
% Check coarseness factor.

legal = abs((3:2:27) - Gfactor) < 4*eps;
if (~any(legal))
  error([' FINE2COARSE: illegal coarseness factor, Gfactor = ',         ...
         num2str(Gfactor)]);
end

% If applicable, get larger grid structure.

if (~isstruct(Ginp))
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

% If extraction indices are not provided, set (Imin,Jmin) and (Imax,Jmax)
% to the maximum coarseness grid possible at PSI-points, such that
% averaging from fine-to-coarse is bounded.

half = (Gfactor-1)/2;

if (numel(varargin) == 0)
  Imin = 1+half;                 % mininum I-value possible (PSI-points)
  Istr = Imin+1:Gfactor:L;
  Iend = Istr+Gfactor-1;  
  Indx = find(Iend > L-half);
  if (~isempty(Indx))
    Istr(Indx)=[];
    Iend(Indx)=[];
  end
  Imax = Iend(end);              % maximum I-value possible (PSI-points)

  Jmin = 1+half;                 % mininum J-value possible (PSI-points)
  Jstr = Jmin+1:Gfactor:M;
  Jend = Jstr+Gfactor-1;  
  Jndx = find(Jend > M-half);
  if (~isempty(Jndx))
    Jstr(Jndx)=[];
    Jend(Jndx)=[];
  end
  Jmax = Jend(end);              % maximum J-value possible (PSI-points)
end

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

if (numel(varargin) > 0)
  Istr = Imin+1:Gfactor:L;
  Iend = Istr+Gfactor-1;  
  Indx = find(Iend > L-half);
  if (~isempty(Indx))
    Istr(Indx)=[];
    Iend(Indx)=[];
  end

  Jstr = Jmin+1:Gfactor:M;
  Jend = Jstr+Gfactor-1;  
  Jndx = find(Jend > M-half);
  if (~isempty(Jndx))
    Jstr(Jndx)=[];
    Jend(Jndx)=[];
  end
end

disp(blanks(1));
disp('Extracting Coarse Grid from Finer Grid: ');
disp(blanks(1));
disp(['  Gfactor = ', num2str(Gfactor)]);
disp(['     Imin = ', num2str(Imin)]);
disp(['     Imax = ', num2str(Imax)]);
disp(['     Jmin = ', num2str(Jmin)]);
disp(['     Jmax = ', num2str(Jmax)]);
disp(blanks(1));
disp('Legal Imax values given Imin:');
disp(blanks(1));
string = sprintf('%5i %5i %5i %5i %5i %5i %5i %5i %5i %5i\n', Iend);
disp(string);
disp(blanks(1));
disp('Legal Jmax values given Jmin:');
disp(blanks(1));
string = sprintf('%5i %5i %5i %5i %5i %5i %5i %5i %5i %5i\n', Jend);
disp(string);

% Now that all the legal values have been reported, remove those values
% exceeding Imax and Jmax.

Indx = find(Iend > Imax);
if (~isempty(Indx))
  Istr(Indx)=[];
  Iend(Indx)=[];
end

Jndx = find(Jend > Jmax);
if (~isempty(Jndx))
  Jstr(Jndx)=[];
  Jend(Jndx)=[];
end

% Check indices to process.

if (numel(varargin) > 0)

  if (Imin >= Imax)
    error([' FINE2COARSE: Imin >= Imax,    ',                           ...
           ' Imin = ', num2str(Imin), ',',                              ...
           ' Imax = ', num2str(Imax)]);
  end

  if (Jmin >= Jmax)
    error([' FINE2COARSE: Jmin >= Jmax,    ',                           ...
           ' Jmin = ', num2str(Jmin), ','                               ...
           ' Jmax = ', num2str(Jmax)]);
  end

  if (Imin < half+1)
    error([' FINE2COARSE: illegal Imin = ', num2str(Imin),              ...
           ',   minimum value possible = ', num2str(half+1)]);
  end

  if (~any(Iend == Imax))
    error([' FINE2COARSE: illegal Imax = ', num2str(Imax),              ...
           ',   maximum value possible = ', num2str(Iend(end)),         ...
           ',   for given Imin = ', num2str(Imin)]);
  end

  if (Jmin < half+1)
    error([' FINE2COARSE: illegal Jmin = ', num2str(Jmin),              ...
           ',   minimum value possible = ', num2str(half+1)]);
  end

  if (~any(Jend == Jmax))
    error([' FINE2COARSE: illegal Jmax = ', num2str(Jmax),              ...
           ',   maximum value possible = ', num2str(Jend(end)),         ...
           ',   given Jmin = ', num2str(Jmin)]);
  end
end

% Set grid variables to process.

grd_vars = {'h', 'f', 'angle', 'pm', 'pn'};

if (F.curvilinear)
  grd_vars = [grd_vars, 'dmde', 'dndx'];
end

grd_vars = [grd_vars, 'x_rho', 'y_rho', 'x_psi', 'y_psi',               ...
                      'x_u', 'y_u', 'x_v', 'y_v'];

if (F.spherical)
  grd_vars = [grd_vars, 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',     ...
                        'lon_u', 'lat_u', 'lon_v', 'lat_v'];
end

grd_vars = [grd_vars, 'mask_rho', 'mask_psi', 'mask_u', 'mask_v'];

grd_vars = [grd_vars, 'hraw'];               % needs to be last

%--------------------------------------------------------------------------
% Extract coaser grid.  Only valid for C-grid...
%--------------------------------------------------------------------------

% Set fine and coarse grid ROMS I-indices to extract.

IindexP = false([1 L  ]);
IindexR = false([1 L+1]);
IindexU = false([1 L  ]);
IindexV = false([1 L+1]);

IpC = Imin     :Gfactor:Imax;
IrC = Imin-half:Gfactor:Imax+Gfactor;
IuC = Imin     :Gfactor:Imax;
IvC = Imin-half:Gfactor:Imax+Gfactor;

IindexP(IpC) = true;
IindexR(IrC) = true;
IindexU(IuC) = true;
IindexV(IvC) = true;

% Set fine and coarse grid ROMS J-indices to extract.

JindexP = false([1 M  ]);
JindexR = false([1 M+1]);
JindexU = false([1 M+1]);
JindexV = false([1 M  ]);

JpC = Jmin     :Gfactor:Jmax;
JrC = Jmin-half:Gfactor:Jmax+Gfactor;
JuC = Jmin-half:Gfactor:Jmax+Gfactor;
JvC = Jmin     :Gfactor:Jmax;

JindexP(JpC) = true;
JindexR(JrC) = true;
JindexU(JuC) = true;
JindexV(JvC) = true;

% Set fine and coarse grids (XI,ETA) coordinates.

[YpF, XpF] = meshgrid(1.0:1:M     , 1.0:1:L     );
[YrF, XrF] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);
[YuF, XuF] = meshgrid(0.5:1:Mp-0.5, 1.0:1:L     );
[YvF, XvF] = meshgrid(1.0:1:M     , 0.5:1:Lp-0.5);

XpC = XpF(IindexP,JindexP);   YpC = YpF(IindexP,JindexP);   
XrC = XrF(IindexR,JindexR);   YrC = YrF(IindexR,JindexR);   
XuC = XuF(IindexU,JindexU);   YuC = YuF(IindexU,JindexU);   
XvC = XvF(IindexV,JindexV);   YvC = YvF(IindexV,JindexV);   

% If requested, plot fine and coarse grid points.

if (Lplot)
  figure;
  plot(XrF, YrF, 'bo',  XrC,  YrC,  'r+',                               ...
       XpF, YpF, 'b--', XpF', YpF', 'b--',                              ...
       XpC, YpC, 'r-',  XpC', YpC', 'r-')
  title('RHO-points in (XI,ETA) Coordinates  (PSI-points mesh)');
  xlabel(['Coarseness Factor = ', num2str(Gfactor), blanks(4),          ...
          '(blue circles: fine RHO-points,',                            ...
	  ' red plus: coarse RHO-points)']);

  figure;
  plot(XuF, YuF, 'bo',  XuC,  YuC,  'r+',                               ...
       XvF, YvF, 'b--', XvF', YvF', 'b--',                              ...
       XvC, YvC, 'r-',  XvC', YvC', 'r-')
  title('U-points in (XI,ETA) Coordinates  (V-points mesh)');
  xlabel(['Coarseness Factor = ', num2str(Gfactor), blanks(4),          ...
          '(blue circles: fine U-points,',                              ...
	  ' red plus: coarse RHO-points)']);

  figure;
  plot(XvF, YvF, 'bo',  XvC,  YvC,  'r+',                               ...
       XvF, YvF, 'b--', XvF', YvF', 'b--',                              ...
       XuC, YuC, 'r-',  XuC', YuC', 'r-')
  title('V-points in (XI,ETA) Coordinates  (U-points mesh)');
  xlabel(['Coarseness Factor = ', num2str(Gfactor), blanks(4),          ...
          '(blue circles: fine V-points,',                              ...
	  ' red plus: coarse RHO-points)']);
end

% Initialize several subsomain structure paameters.

C.grid_name = Gout;
C.Lm = length(IpC)-1;
C.Mm = length(JpC)-1;
C.spherical = F.spherical;

if isfield(F, 'uniform')
  C.uniform = F.uniform;
else
  C.uniform = 0;
end

% Extract grid variables.

for value = grd_vars
  field = char(value);
  got.(field) = false;
  if (isfield(F,field))
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

if (got.x_psi && got.y_psi)
  C.xl = max(C.x_psi(:)) - min(C.x_psi(:));
  C.el = max(C.y_psi(:)) - min(C.y_psi(:));
else
  C.xl = 0;
  C.el = 0;
end

% Process inverse grid spacing array "pm" and "pn".  Sum finer grid
% values to conserve area.

Ibry  = Istr(1)-(half+1):Gfactor:Iend(end)+(half+1);
Jbry  = Jstr(1)-(half+1):Gfactor:Jend(end)+(half+1);

dxC = nan(size(XrC));
dyC = nan(size(XrC));

% X-spacing interior points including southern and northern boundaries.
% The values are summed for each finer grid j-value within the coarser
% grid cell and then averaged at the end.

for j = 1:length(Jbry)
  Javg = Jbry(j)-half:1:Jbry(j)+half;
  indx = Javg < 1 | Javg > M;
  Javg(indx) = [];
  Nvalues = length(Javg);
  values  = zeros([1 Nvalues]);
  for i = 1:length(Istr)
    for jc = 1:Nvalues
      values(jc) = sum(1 ./ squeeze(F.pm(Istr(i):Iend(i), Javg(jc))));
    end
    dxC(i+1,j) = mean(values(:));
  end
end

% Y-spacing interior points including western and eastern boundaries.
% The values are summed for each finer grid j-value within the coarser
% grid cell and averaged at the ne

for i = 1:length(Ibry)
  Iavg = Ibry(i)-half:1:Ibry(i)+half;
  indx = Iavg < 1 | Iavg > L;
  Iavg(indx) = [];
  Nvalues = length(Iavg);
  values  = zeros([1 Nvalues]);
  for j = 1:length(Jstr)
    for ic = 1:Nvalues
      values(ic) = sum(1 ./ squeeze(F.pn(Iavg(ic), Jstr(j):Jend(j))));
    end
    dyC(i,j+1) = mean(values(:));
  end
end

%  Boundary points.  The values depend on the extraction indices.  We
%  either have enough points for averaging or we need to mirror values
%  by multypling by 2.

%  Western boundary.

if (Imin < Gfactor)
  Iwest  = Imin-half:1:Imin;
  factor = 0.5;
else
  Iwest = Imin-(Gfactor-1):1:Imin;
  factor = 1;
end  

for j = 1:length(Jbry)
  Javg = Jbry(j)-half:1:Jbry(j)+half;
  indx = Javg < 1 | Javg > M;
  Javg(indx) = [];
  Nvalues = length(Javg);
  values = zeros([1 Nvalues]);
  for jc = 1:Nvalues
    val  = 1 ./ squeeze(F.pm(Iwest, Javg(jc)));
    val(1) = factor * val(1);
    values(jc) = sum(val(:)) / factor;
  end      
  dxC(1,j) = mean(values(:));
end

%  Eastern boundary.

if (Imax > (Lp-Gfactor))
  Ieast = Iend(end)+1:1:Iend(end)+1+half;
  factor = 0.5;
else
  Ieast = Iend(end)+1:1:Iend(end)+Gfactor;
  factor = 1;
end

for j = 1:length(Jbry)
  Javg = Jbry(j)-half:1:Jbry(j)+half;
  indx = Javg < 1 | Javg > M;
  Javg(indx) = [];
  Nvalues = length(Javg);
  values = zeros([1 Nvalues]);
  for jc = 1:Nvalues
    val = 1 ./ squeeze(F.pm(Ieast, Javg(jc)));
    val(end) = factor * val(end);
    values(jc) = sum(val(:)) / factor;
  end
  dxC(end,j) = mean(values(:));
end

%  Southern boundary.

if (Jmin < Gfactor)
  Jsouth = Jmin-half:1:Jmin;
  factor = 0.5;
else
  Jsouth = Jmin-(Gfactor-1):1:Jmin;
  factor = 1;
end

for i = 1:length(Ibry)
  Iavg = Ibry(i)-half:1:Ibry(i)+half;
  indx = Iavg < 1 | Iavg > L;
  Iavg(indx) = [];
  Nvalues = length(Iavg);
  values  = zeros([1 Nvalues]);
  for ic = 1:Nvalues  
    val = 1 ./ squeeze(F.pn(Iavg(ic), Jsouth));
    val(1) = factor * val(1);
    values(ic) = sum(val(:)) / factor;
  end
  dyC(i,1) = mean(values(:));
end

%  Northern boundary.

if (Jmax > (Mp-Gfactor))
  Jnorth = Jend(end)+1:1:Jend(end)+1+half;
  factor = 0.5;
else
  Jnorth = Jend(end)+1:1:Jend(end)+Gfactor;
  factor = 1;
end  
  
for i = 1:length(Ibry)
  Iavg = Ibry(i)-half:1:Ibry(i)+half;
  indx = Iavg < 1 | Iavg > L;
  Iavg(indx) = [];
  Nvalues = length(Iavg);
  values  = zeros([1 Nvalues]);
  for ic = 1:Nvalues
    val = 1 ./ squeeze(F.pn(Iavg(ic), Jnorth));
    val(end) = factor * val(end);
    values(ic) = sum(val(:)) / factor;
  end
  dyC(i,end) = mean(values(:));
end

C.pm = 1 ./ dxC;
C.pn = 1 ./ dyC;

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

if (got.hraw)
  bath = size(C.hraw,3);
  for rec=1:bath
    status = nc_write (Gout, 'hraw', C.hraw, rec);
    if (status ~= 0), return, end
  end
else
  status = nc_write (Gout, 'hraw', C.h, 1);
  if (status ~= 0), return, end
end

for value = 1:length(grd_vars)-1
  field = char(grd_vars(value));
  if (got.(field))
    status = nc_write (Gout, field, C.(field));
    if (status ~= 0), return, end
  end  
end

%--------------------------------------------------------------------------
% Add coastline data if available.
%--------------------------------------------------------------------------

if (isfield(F, 'lon_coast') && isfield(F, 'lat_coast'))
  add_coastline (Gout, F.lon_coast, F.lat_coast);
end

%--------------------------------------------------------------------------
% Get full extracted grid structure.
%--------------------------------------------------------------------------

C = get_roms_grid(Gout);

return
