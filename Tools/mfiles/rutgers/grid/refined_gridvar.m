function [Fvalues,S] = refined_gridvar(Cinp, Finp, field, varargin)

%
% REFINED_GRIDVAR: Sets requested refined grid field from ROMS coarse grid
%
% [Fvalues,S] = refined_gridvar(Cinp, Finp, field, DivideCoarse)
%
% Sets requested refined grid variable from ROMS coarse grid division or
% quadratic conservative interpolation (Clark and Farley, 1984), such
% that the coarse-to-fine and fine-to-coarse are reversible:
%
%       AVG[F.field(:,:)_i , i=1,rfactor] = C.field(Idg,Jdg)
%
% where (Idg,Jdg) are the indices of the coarse grid donor cell contatining
% the finer grid.
%
% If field to process is uniform, the interpolation is not performed to
% avoid round-off. It should be only used to process 2D ROMS grid fields
% like: 'angle', 'dmde', 'dmdx', 'f', 'h', 'pm', and 'pn'.
%
% On Input:
%
%    Cinp          ROMS coarse Grid NetCDF files/URL name (string)
%              or, an existing ROMS coarse grid structure (struct array)
%
%    Finp          ROMS finer Grid NetCDF files/URL name (string)
%              or, an existing ROMS finer grid structure (struct array)
%
%    field         Field name to process (NetCDF variable name or
%                  structure field)
%
%    DivideCoarse  Switch to return the coarse value divided by the grid
%                  refinement ratio.  This only apply for 'pm' or 'pn'
%                  fields to garanttee exact area conservation (OPTIONAL,
%                  default = false)
%
% On Output:
%
%    Fvalues       Interpolated field (2D Array)
%
%    S             Interpolation information structure (struct, OPTIONAL):
%
%                    S.grid_time       Arakawa grid C-type (PSI, RHO, u, v)
%                    S.is2d            processing 2D field (false for 3D)
%                    S.spherical       spherical switch
%                    S.uniform         uniform field switch
%                    S.masked          land masked weights switch
%                    S.refine_factor   grid refinement factor
%                    S.parent_Imin     coarse grid lower-left  I-coordinate
%                    S.parent_Imax     coarse grid upper-right I-coordinate
%                    S.parent_Jmin     coarse grid lower-left  J-coordinate
%                    S.parent_Jmax     coarse grid upper-right J-coordinate
%                    S.x               Cartesian or spherical X-coordinate
%                    S.y               Cartesian or spherical Y-coordinate
%                    S.val             interpolated values
%                    S.p               interpolation p-fractional distance
%                    S.q               interpolation q-fractional distance
%                    S.w               quadratic interpolation weights
%                    S.xavg            averaged field interior X-coordinate
%                    S.yavg            averaged field interior Y-coordinate
%                    S.avg             averaged field values: fine2coarse
%                    S.sum             accumulate field values: fine2coarse
%
% You can plot interpolation fields as follows:
%
%  pcolor(S.x, S.y, V); shading facetted; colorbar;
%  pcolor(S.x, S.y, S.val); shading facetted; colorbar;   % same as above
%
%  pcolor(S.xavg, S.yavg, S.avg); shading facetted; colorbar;
%  pcolor(S.xavg, S.yavg, S.sum); shading facetted; colorbar;
%
  
% svn $Id: refined_gridvar.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.
  
masked = false;         % interpolation weights are not scaled by land mask
isgridspacing = false;  % processing 'pm' or 'pn'
got_xy = false;         % Cartesian coodinates (x_rho, y_rho) ...

% Set optional parameters.

DivideCoarse = false;

switch numel(varargin)
  case 1
    DivideCoarse = varargin{1};
end

%--------------------------------------------------------------------------
% Set coarse and fine grids structures.
%--------------------------------------------------------------------------

if (~isstruct(Cinp))
  C = get_roms_grid(Cinp);
else
  C = Cinp;
end

if (~isstruct(Finp))
  F = get_roms_grid(Finp);
else
  F = Finp;
end

if (~isfield(C,field))
  error(['REFINED_GRIDVAR: cannot find field in coarse grid structure,' ...
         ' field = ', field]);
end

% Get finer grid extraction indices (PSI-points) and refinement factor.

if (isfield(F,'parent_Imin'))
  Imin = F.parent_Imin;
else
  error(['REFINED_GRIDVAR: cannot find field in finer grid structure:', ...
         ' parent_Imin']);
end

if (isfield(F,'parent_Imax'))
  Imax = F.parent_Imax;
else
  error(['REFINED_GRIDVAR: cannot find field in finer grid structure:', ...
         ' parent_Imax']);
end

if (isfield(F,'parent_Jmin'))
  Jmin = F.parent_Jmin;
else
  error(['REFINED_GRIDVAR: cannot find field in finer grid structure:', ...
         ' parent_Jmin']);
end

if (isfield(F,'parent_Jmax'))
  Jmax = F.parent_Jmax;
else
  error(['REFINED_GRIDVAR: cannot find field in finer grid structure:', ...
         ' parent_Jmax']);
end

if (isfield(F,'refine_factor'))
  rfactor = F.refine_factor;
else
  error(['REFINED_GRIDVAR: cannot find field in finer grid structure:', ...
         ' refine_factor']);
end

% Set various flags.

switch field
  case {'pm', 'pn'}
    Ctype = 'rho_point';
    isgridspacing = true;
    if (all(isfield(F,{'x_rho', 'y_rho', 'x_u', 'y_u', 'x_v', 'y_v'})))
      if ~(isempty(F.x_rho) && isempty(F.y_rho) &&                      ...
           isempty(F.x_u  ) && isempty(F.y_u  ) &&                      ...
           isempty(F.x_v  ) && isempty(F.y_v  ))
        got_xy = true;
      end
    end
  case {'angle', 'dmde', 'dndx', 'f', 'h', 'hraw'}
    Ctype = 'rho_point';
  otherwise
    error([' REFINED_GRIDVAR: illegal field: ', field])
end

% Determine if the interpolated value desired is just the coarse grid
% value divided grid refined factor.  This is only appliable to the
% inverse grid spacing variables 'pm' and 'pn'.  This can be done to
% guaranttee that the coarse grid area is mantained in the finer grid.

DivideCoarseValue = DivideCoarse & isgridspacing;

% Check if field to process is uniform.

isuniform = (length(unique(C.(field))) == 1);

%--------------------------------------------------------------------------
% Interpolate requested field.
%--------------------------------------------------------------------------

% Set coordinates.

switch Ctype
  case 'rho_point'
    if (C.spherical)
      XF = F.lon_rho;
      YF = F.lat_rho;
    else
      XF = F.x_rho;
      YF = F.y_rho;
    end
  case 'u_point'
    if (C.spherical)
      XF = F.lon_u;
      YF = F.lat_u;
    else
      XF = F.x_u;
      YF = F.y_u;
    end
  case 'v_point'
    if (C.spherical)
      XF = F.lon_v;
      YF = F.lat_v;
    else
      XF = F.x_v;
      YF = F.y_v;
    end
end

[Lp,Mp] = size(XF);

L = Lp-1;   Lm = L-1;
M = Mp-1;   Mm = M-1;

% If uniform field, interpolation is not needed.  This is done to avoid
% round-off due to high order interpolation.  Notice that for 'pm' and
% 'pn' we need to multipy by rfactor since dx=1/pm and dy=1/pn.

if (isuniform)

% Interpolating an uniform (constant) field.

  switch field
    case {'pm', 'pn'}
      constant_val = unique(C.(field)) * rfactor;
    otherwise
      constant_val = unique(C.(field));
  end

  [Im,Jm] = size(XF);
  Fvalues = ones([Im Jm]) .* constant_val;

else

% Interpolating a non uniform field.

  delta = 1.0/rfactor;
  half  = 0.5*delta;

  IpF = (Imin:delta:Imax);                           % PSI-points
  JpF = (Jmin:delta:Jmax);                           % PSI-points
  IrF = [IpF(1)-half IpF+half];                      % RHO-points
  JrF = [JpF(1)-half JpF+half];                      % RHO-points

  switch Ctype
    case 'rho_point'
      switch field
        case {'pm'}
          if (got_xy)
            dx = zeros(size(XF));
            dx(2:L,1:Mp) = sqrt((F.x_u  (2:L ,1:Mp) -                   ...
                                 F.x_u  (1:Lm,1:Mp)).^2 +               ...
                                (F.y_u  (2:L ,1:Mp) -                   ...
                                 F.y_u  (1:Lm,1:Mp)).^2);
            dx(1  ,1:Mp) = sqrt((F.x_u  (1   ,1:Mp) -                   ...
                                 F.x_rho(1   ,1:Mp)).^2 +               ...
                                (F.y_u  (1   ,1:Mp) -                   ...
                                 F.y_rho(1   ,1:Mp)).^2).*2.0;
            dx(Lp ,1:Mp) = sqrt((F.x_rho(Lp  ,1:Mp) -                   ...
                                 F.x_u  (L   ,1:Mp)).^2 +               ...
                                (F.y_rho(Lp  ,1:Mp) -                   ...
                                 F.y_u  (L   ,1:Mp)).^2).*2.0;
            Fvalues = 1 ./ dx;
          else
            Fvalues = C.(field)(fix(IrF),fix(JrF)) .* rfactor;
          end       
        case {'pn'}
          if (got_xy)
            dy = zeros(size(XF));
            dy(1:Lp,2:M) = sqrt((F.x_v  (1:Lp,2:M ) -                   ...
                                 F.x_v  (1:Lp,1:Mm)).^2 +               ...
                                (F.y_v  (1:Lp,2:M ) -                   ...
                                 F.y_v  (1:Lp,1:Mm)).^2);
            dy(1:Lp,1  ) = sqrt((F.x_v  (1:Lp,1   ) -                   ...
                                 F.x_rho(1:Lp,1   )).^2 +               ...
                                (F.y_v  (1:Lp,1   ) -                   ...
                                 F.y_rho(1:Lp,1   )).^2).*2.0;
            dy(1:Lp,Mp ) = sqrt((F.x_rho(1:Lp,Mp  ) -                   ...
                                 F.x_v  (1:Lp,M   )).^2 +               ...
                                (F.y_rho(1:Lp,Mp  ) -                   ...
                                 F.y_v  (1:Lp,M   )).^2).*2.0;          
            Fvalues = 1 ./ dy;
          else
            Fvalues = C.(field)(fix(IrF),fix(JrF)) .* rfactor;
          end
        otherwise
          [Fvalues, S] = qc_interp(C, F, C.(field), masked);
      end
  end
  [Im,Jm] = size(XF);

end

%--------------------------------------------------------------------------
% Set optional output arguments.
%--------------------------------------------------------------------------

if (nargout > 1)
  if ~(isuniform || DivideCoarseValue)  
    W.field = field;
    names = [fieldnames(W); fieldnames(S)];
    S = cell2struct([struct2cell(W); struct2cell(S)], names, 1);
  else
    S.field = field;
    S.grid_type = Ctype;
    S.is2d = true;
    S.spherical = C.spherical;
    S.isuniform = isuniform;
    S.masked = masked;
    S.refine_factor = rfactor;
    S.parent_Imin = Imin;
    S.parent_Imax = Imax;
    S.parent_Jmin = Jmin;
    S.parent_Jmax = Jmax;
    S.x = XF;
    S.y = YF;
    S.val = Fvalues;
    S.p = [];
    S.q = [];
    S.w = [];

% Average interpolated field interior points back to coarse grid locations.
% This equivalent to the fine-to-coarse operation.

    Imid = 2+(rfactor-1)/2:rfactor:Im;
    Jmid = 2+(rfactor-1)/2:rfactor:Jm;

    S.xavg = S.x(Imid, Jmid);
    S.yavg = S.y(Imid, Jmid);

    Istr = 2:rfactor:Im-1;
    Iend = Istr(2)-1:rfactor:Im;

    Jstr = 2:rfactor:Jm-1;
    Jend = Jstr(2)-1:rfactor:Jm;

    S.avg = nan([length(Imid) length(Jmid)]);
    S.sum = nan(size(S.avg));
  
    for j = 1:length(Jstr)
      for i = 1:length(Istr)
        values = S.val(Istr(i):Iend(i), Jstr(j):Jend(j));
        S.avg(i,j) = mean(values(:));
        S.sum(i,j) = sum(values(:));
      end
    end
  end
end

return

