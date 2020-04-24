function [Fvalues,S]=qc_interp(Cinp, Finp, Cvalues, Lmask)

%
% QC_INTERP: Quadratic conservative interpolation for ROMS refinement
%
% [Fvalues,S]=qc_interp(Cinp, Finp, Cvalues, Lmask)
%
% Interpolates requested field from ROMS coarse to fine grid using
% 9-points quadratic conservative interpolation (Clark and Farley, 1984),
% such that the coarse-to-fine and fine-to-coarse are reversible:
%
%       AVG[Fvalues(:,:)_i , i=1,rfactor^2] = Cvalues(Id,Jd)
%
% where (Id,Jd) are the indices of the coarse grid donor cell contatining
% the finer grid values (see diagram below).
%
% This function can also be use to interpolate ROMS 3D initial conditions.
% Currently, no vertical interpolation is carried out since in refinement
% the coarse and finer grid has the same number of vertical levels.  This
% constraint may be relaxed in the future.
%
% On Input:
%
%    Cinp          ROMS coarse Grid NetCDF files/URL name (string)
%              or, an existing ROMS coarse grid structure computed with
%                  'get_roms_grid' (struct)
%
%    Finp          ROMS finer Grid NetCDF files/URL name (string)
%              or, an existing ROMS finer grid structure computed with
%                  'get_roms_grid' (struct)
%
%    Cvalues       Coarse field values to interpolated from (2D or 3D
%                  array at any C-grid location)
%
%    Lmask         Switch to adjust the quadratic interpolation weights
%                  by the land mask.
%
% On Output:
%
%    Fvalues       Interpolated field values (2D or 3D Array)
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
%
% The following diagram show the conventions for quadratic conservative
% interpolation for a receiver grid with coordinates (Ir,Jr) from the donor
% cells (Id-1, Id, Id+1; Jd-1, Jd, Jd+1):
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
% where Rm, Ro, Rp and Sm, So, Sp are the interpolant functions in the
% I- and J-directions, respectively:
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
%                                                   (no refinement)
% The quadratic interpolation is:
%
%     Fvalues(Ir, Jr) = Rm * Sm * C(Id-1, Jd-1) +
%                       Ro * Sm * C(Id  , Jd-1) +
%                       Rp * Sm * C(Id+1, Jd-1) +
%                       Rm * So * C(Id-1, Jd  ) +
%                       Ro * So * C(Id  , Jd  ) +
%                       Rp * So * C(Id+1, Jd  ) +
%                       Rm * Sp * C(Id-1, Jd+1) +
%                       Ro * Sp * C(Id  , Jd+1) +
%                       Rp * Sp * C(Id+1, Jd+1)
%

% svn $Id: qc_interp.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set coarse and fine grid structures.

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

% Get finer grid extraction indices (PSI-points) and refinement factor.

if (isfield(F,'parent_Imin'))
  Imin = F.parent_Imin;
else
  error([' QC_INTERP: cannot find field in finer grid structure:',      ...
         ' parent_Imax']);
end

if (isfield(F,'parent_Imax'))
  Imax = F.parent_Imax;
else
  error([' QC_INTERP: cannot find field in finer grid structure:',      ...
         ' parent_Imax']);
end

if (isfield(F,'parent_Jmin'))
  Jmin = F.parent_Jmin;
else
  error([' QC_INTERP: cannot find field in finer grid structure:',      ...
         ' parent_Jmin']);
end

if (isfield(F,'parent_Jmax'))
  Jmax = F.parent_Jmax;
else
  error([' QC_INTERP: cannot find field in finer grid structure:',      ...
         ' parent_Jmax']);
end

if (isfield(F,'refine_factor'))
  rfactor = F.refine_factor;
else
  error([' QC_INTERP: cannot find field in finer grid structure:',      ...
         ' refine_factor']);
end

% Set coarse grid size.

L = C.Lm+1;  Lp = L+1;
M = C.Mm+1;  Mp = M+1;

% Check input field grid type.

Crank = length(size(Cvalues));
is2d = false;
is3d = false;

if (Crank == 2)
  is2d = true;
  [Im,Jm] = size(Cvalues);
elseif (Crank == 3)
  is3d = true;
  [Im,Jm,Km] = size(Cvalues);
else
  error([' QC_INTERP: input field must be a 2D or 3D array, rank = ',   ...
         num2str(Crank)]);
end

if (Im == L && Jm == M)
  Ctype = 'psi_point';
elseif (Im == Lp && Jm == Mp)
  Ctype = 'rho_point';
elseif (Im == L && Jm == Mp)
  Ctype = 'u_point';
elseif (Im == Lp && Jm == M)
  Ctype = 'v_point';
else
  error([' QC_INTERP: unable to determine Arakawa C-grid type:',        ...
         ' Cvalues is inconsistent with coarse grid structure']);
end

% Initialize internal switch for masking interpolation weights.

masked = false;

%--------------------------------------------------------------------------
% Interpolate to finer grid.
%--------------------------------------------------------------------------

% Set finer grid fractional coordinates.

delta = 1.0/rfactor;
half  = 0.5*delta;

IpF = (Imin:delta:Imax);                             % PSI-points
JpF = (Jmin:delta:Jmax);                             % PSI-points
IrF = [IpF(1)-half IpF+half];                        % RHO-points
JrF = [JpF(1)-half JpF+half];                        % RHO-points

% Get linear indices of the coarse grid cells containing the finer grid
% cells.  The I and J indices are shifted by one since Matlab does not
% support zero index arrays.

switch Ctype
  case 'psi_point'
    [YC, XC] = meshgrid(1.0:1:M     , 1.0:1:L     ); 
    [YF, XF] = meshgrid(JpF, IpF);
    [Im, Jm] = size(XC);

    Idg   = fix(XF - 0.5);                                % Idg
    Idgm1 = max(1, min(Iu, Idg-1));                       % Idg-1
    Idgp1 = min(Iu, Idg+1);                               % Idg+1

    Jdg   = fix(YF - 0.5);                                % Jdg
    Jdgm1 = max(1, min(Jm, Jdg-1));                       % Jdg-1
    Jdgp1 = min(Jm, Jdg+1);                               % Jdg+1

    if (Lmask)
      if (isfield(C,'mask_psi'))
        mask = C.mask_psi;
      else
        mask = ones([Im Jm]);
      end      
    end
  case 'rho_point'
    [YC, XC] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);
    [YF, XF] = meshgrid(JrF, IrF);
    [Im, Jm] = size(XC);

    Idg   = fix(XF - 0.5) + 1;                            % Idg
    Idgm1 = max(1, min(Im, Idg-1));                       % Idg-1
    Idgp1 = min(Im, Idg+1);                               % Idg+1

    Jdg   = fix(YF - 0.5) + 1;                            % Jdg
    Jdgm1 = max(1, min(Jm, Jdg-1));                       % Jdg-1
    Jdgp1 = min(Jm, Jdg+1);                               % Jdg+1

    if (Lmask)
      if (isfield(C,'mask_rho'))
        mask = C.mask_rho;
      else
        mask = ones([Im Jm]);
      end      
    end
  case 'u_point'
    [YC, XC] = meshgrid(0.5:1:Mp-0.5, 1.0:1:L     ); 
    [YF, XF] = meshgrid(JrF, IpF);
    [Im, Jm] = size(XC);

    Idg   = fix(XF - 0.5);                                % Idg
    Idgm1 = max(1, min(Iu, Idg-1));                       % Idg-1
    Idgp1 = min(Iu, Idg+1);                               % Idg+1

    Jdg   = fix(YF - 0.5) + 1;                            % Jdg
    Jdgm1 = max(1, min(Ju, Jdg-1));                       % Jdg-1
    Jdgp1 = min(Ju, Jdg+1);                               % Jdg+1

    if (Lmask)
      if (isfield(C,'mask_u'))
        mask = C.mask_u;
      else
        mask = ones([Im Jm]);
      end      
    end
  case 'v_point'
    [YC, XC] = meshgrid(1.0:1:M     , 0.5:1:Lp-0.5);
    [YF, XF] = meshgrid(JpF, IrF);
    [Im, Jm] = size(XC);

    Idg   = fix(XF - 0.5) + 1;                            % Idg
    Idgm1 = max(1, min(Im, Idg-1));                       % Idg-1
    Idgp1 = min(Im, Idg+1);                               % Idg+1
    
    Jdg   = fix(YF - 0.5);                                % Jdg
    Jdgm1 = max(1, min(Jm, Jdg-1));                       % Jdg-1
    Jdgp1 = min(Jm, Jdg+1);                               % Jdg+1

    if (Lmask)
      if (isfield(C,'mask_v'))
        mask = C.mask_v;
      else
        mask = ones([Im Jm]);
      end      
    end
  otherwise
    [YC, XC] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);
    [YF, XF] = meshgrid(JrF, IrF);
    [Im, Jm] = size(XC);

    Idg   = fix(XF - 0.5) + 1;                            % Idg
    Idgm1 = max(1, min(Im, Idg-1));                       % Idg-1
    Idgp1 = min(Im, Idg+1);                               % Idg+1

    Jdg   = fix(YF - 0.5) + 1;                            % Jdg
    Jdgm1 = max(1, min(Jm, Jdg-1));                       % Jdg-1
    Jdgp1 = min(Jm, Jdg+1);                               % Jdg+1

    if (Lmask)
      if (isfield(C,'mask_rho'))
        mask = C.mask_rho;
      else
        mask = ones([Im Jm]);
      end      
    end
end

% Set coarse grid indices containing refinement points.

index1 = sub2ind([Im,Jm], Idgm1, Jdgm1);                  % (Idg-1, Jdg-1)
index2 = sub2ind([Im,Jm], Idg  , Jdgm1);                  % (Idg  , Jdg-1)
index3 = sub2ind([Im,Jm], Idgp1, Jdgm1);                  % (Idg+1, Jdg-1)
index4 = sub2ind([Im,Jm], Idgm1, Jdg  );                  % (Idg-1, Jdg  )
index5 = sub2ind([Im,Jm], Idg  , Jdg  );                  % (Idg  , Jdg  )
index6 = sub2ind([Im,Jm], Idgp1, Jdg  );                  % (Idg+1, Jdg  )
index7 = sub2ind([Im,Jm], Idgm1, Jdgp1);                  % (Idg-1, Jdg+1)
index8 = sub2ind([Im,Jm], Idg  , Jdgp1);                  % (Idg  , Jdg+1)
index9 = sub2ind([Im,Jm], Idgp1, Jdgp1);                  % (Idg+1, Jdg+1)    

% Set fractional distances.

p = (XF - XC(index5)) ./ (XC(index6) - XC(index5));
q = (YF - YC(index5)) ./ (YC(index8) - YC(index5));

if (any(isnan(p)))
  p(isnan(p)) = 0;
end

if (any(isnan(q)))
  q(isnan(q)) = 0;
end

% Compute quadratic interpolation weights.

if (rfactor > 0)
  alpha = ((1 / rfactor)^2 - 1) / 24;
else
  alpha = 0;
end

Rm = 0.5 .* p .* (p - 1) + alpha;
Ro = (1 - p .* p)        - 2 * alpha;
Rp = 0.5 .* p .* (p + 1) + alpha;

Sm = 0.5 .* q .* (q - 1) + alpha;
So = (1 - q .* q)        - 2 * alpha;
Sp = 0.5 .* q .* (q + 1) + alpha;

if (Lmask)
  w(1,:,:) = mask(index1) .* Rm .* Sm;
  w(2,:,:) = mask(index2) .* Ro .* Sm;
  w(3,:,:) = mask(index3) .* Rp .* Sm;
  w(4,:,:) = mask(index4) .* Rm .* So;
  w(5,:,:) = mask(index5) .* Ro .* So;
  w(6,:,:) = mask(index6) .* Rp .* So;
  w(7,:,:) = mask(index7) .* Rm .* Sp;
  w(8,:,:) = mask(index8) .* Ro .* Sp;
  w(9,:,:) = mask(index9) .* Rp .* Sp;

  wsum = squeeze(sum(w, 1));
  MaskSum = mask(index1) + mask(index2) +  mask(index3) +               ...
            mask(index4) + mask(index5) +  mask(index6) +               ...
            mask(index7) + mask(index8) +  mask(index9);

  ind = find(MaskSum < 9 & wsum > 0);
  if (~isempty(ind))
    masked = true;
    for n=1:9
      w(n,ind) = squeeze(w(n,ind)) ./ wsum(ind);
    end
  end
else
  w(1,:,:) = Rm .* Sm;
  w(2,:,:) = Ro .* Sm;
  w(3,:,:) = Rp .* Sm;
  w(4,:,:) = Rm .* So;
  w(5,:,:) = Ro .* So;
  w(6,:,:) = Rp .* So;
  w(7,:,:) = Rm .* Sp;
  w(8,:,:) = Ro .* Sp;
  w(9,:,:) = Rp .* Sp;
end

% Interpolate requested field.

if (is2d)

  isuniform = (length(unique(Cvalues)) == 1);
  if (isuniform)
    constant_val = unique(Cvalues);
    Fvalues = ones(size(XF)) .* constant_val;
  else
    Fvalues = squeeze(w(1,:,:)) .* Cvalues(index1) +                    ...
              squeeze(w(2,:,:)) .* Cvalues(index2) +                    ...
              squeeze(w(3,:,:)) .* Cvalues(index3) +                    ...
              squeeze(w(4,:,:)) .* Cvalues(index4) +                    ...
              squeeze(w(5,:,:)) .* Cvalues(index5) +                    ...
              squeeze(w(6,:,:)) .* Cvalues(index6) +                    ...
              squeeze(w(7,:,:)) .* Cvalues(index7) +                    ...
              squeeze(w(8,:,:)) .* Cvalues(index8) +                    ...
              squeeze(w(9,:,:)) .* Cvalues(index9);
  end
  
elseif (is3d)

  Fvalues = nan([size(XF) Km]);
  for k=1:Km
    isuniform = (length(unique(squeeze(Cvalues(k,:,:)))) == 1);
    if (isuniform)
      constant_val = unique(squeeze(Cvalues(k,:,:)));
      Fvalues(k,:,:) = ones(size(XF)) .* constant_val;
    else
      Fvalues(k,:,:) = squeeze(w(1,:,:)) .* squeeze(Cvalues(k,index1))+ ...
                       squeeze(w(2,:,:)) .* squeeze(Cvalues(k,index2))+ ...
                       squeeze(w(3,:,:)) .* squeeze(Cvalues(k,index3))+ ...
                       squeeze(w(4,:,:)) .* squeeze(Cvalues(k,index4))+ ...
                       squeeze(w(5,:,:)) .* squeeze(Cvalues(k,index5))+ ...
                       squeeze(w(6,:,:)) .* squeeze(Cvalues(k,index6))+ ...
                       squeeze(w(7,:,:)) .* squeeze(Cvalues(k,index7))+ ...
                       squeeze(w(8,:,:)) .* squeeze(Cvalues(k,index8))+ ...
                       squeeze(w(9,:,:)) .* squeeze(Cvalues(k,index9));
    end
  end

end

%--------------------------------------------------------------------------
% Set optional output arguments.
%--------------------------------------------------------------------------

if (nargout > 1)
  S.grid_type = Ctype;
  S.is2d = is2d;
  S.spherical = C.spherical;
  S.isuniform = isuniform;
  S.masked = masked;
  S.refine_factor = rfactor;
  S.parent_Imin = Imin;
  S.parent_Imax = Imax;
  S.parent_Jmin = Jmin;
  S.parent_Jmax = Jmax;
  switch Ctype
    case 'psi_point'
      if (C.spherical)
        S.x = F.lon_psi;
        S.y = F.lat_psi;
      else
        S.x = F.x_psi;
        S.y = F.y_psi;
      end
    case 'rho_point'
      if (C.spherical)
        S.x = F.lon_rho;
        S.y = F.lat_rho;
      else
        S.x = F.x_rho;
        S.y = F.y_rho;
      end
    case 'u_point'
      if (C.spherical)
        S.x = F.lon_u;
        S.y = F.lat_u;
      else
        S.x = F.x_u;
        S.y = F.y_u;
      end
    case 'v_point'
      if (C.spherical)
        S.x = F.lon_v;
        S.y = F.lat_v;
      else
        S.x = F.x_v;
        S.y = F.y_v;
      end
  end
  S.val = Fvalues;
  S.p = p;
  S.q = q;
  S.w = w;

% Average interpolated field interior points back to coarse grid locations.
% This equivalent to the fine-to-coarse operation.

  [Im,Jm] = size(S.x);
  
  Imid = 2+(rfactor-1)/2:rfactor:Im;
  Jmid = 2+(rfactor-1)/2:rfactor:Jm;

  S.xavg = S.x(Imid, Jmid);
  S.yavg = S.y(Imid, Jmid);

  Istr = 2:rfactor:Im-1;
  Iend = Istr(2)-1:rfactor:Im;

  Jstr = 2:rfactor:Jm-1;
  Jend = Jstr(2)-1:rfactor:Jm;

  if (is2d)
    S.avg = nan([length(Imid) length(Jmid)]);
    S.sum = nan(size(S.avg));
    for j = 1:length(Jstr)
      for i = 1:length(Istr)
        values = S.val(Istr(i):Iend(i), Jstr(j):Jend(j));
        S.avg(i,j) = mean(values(:));
        S.sum(i,j) = sum(values(:));
      end
    end
  elseif (is3d)
    S.avg = nan([length(Imid) length(Jmid) Km]);
    S.sum = nan(size(S.avg));
    for k = 1:Km
      for j = 1:length(Jstr)
        for i = 1:length(Istr)
          values = squeeze(S.val(Istr(i):Iend(i), Jstr(j):Jend(j), k));
          S.avg(i,j,k) = mean(values(:));
          S.sum(i,j,k) = sum(values(:));
        end
      end
    end
  end  
end

return
