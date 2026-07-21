function F = add_uv_agrid(G, ncfile, u_varname, v_varname)

%
% ADD_UV_AGRID:  Adds rotated vector components to ROMS NetCDF file
%
% add_uv_agrid(ncfile, u_varname, v_varname)
%
% Adds vector components at RHO-points (A-grid) and rotate its values
% to geographical Eastward and Nortward directions.  
%
% On Input:
%
%    G             ROMS grid NetCDF filename (string)
%              or, an existing ROMS grid structure (struct array)
%    ncfile        ROMS NetCDF filename (string)
%    u_varname     U-component variable name to add (string)
%    v_varname     V-component variable name to add (string)
%
% Example: Compute and sdd 4D-Var standard deviations at RHO-points
%          (A-grid) rotated to Eastward and Northward directions.
%
%    add_uv_agrid('roms_grd.nc', 'roms_std.nc', 'u_eastward', 'v_northward)
  
% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Initialize.
  
if (~isstruct(G))
  ncgrid = G;
  G = get_roms_grid(ncgrid);
end
  
%--------------------------------------------------------------------------
% Inquire input NetCDF file.
%--------------------------------------------------------------------------

I = nc_inq(ncfile);

got.Uvar = any(strcmp({I.Variables.Name}, u_varname));
got.Vvar = any(strcmp({I.Variables.Name}, v_varname));

if (any(strcmp({I.Variables.Name}, 'spherical')))
  spherical = nc_read(ncfile, 'spherical');
  if (ischar(spherical))
    if (spherical == 'T' || spherical == 't')
      spherical = true;
    else
      spherical = false;
    end
  end
else
  spherical = true;
end

%--------------------------------------------------------------------------
%  If appropriate, define sponge variables.
%--------------------------------------------------------------------------

ic = 0;

append_vars = false;

if (~got.Uvar)
  ic = ic + 1;
  S.Dimensions = I.Dimensions;
  S.Variables(ic) = roms_metadata(u_varname, spherical);
  append_vars = true;
else
  disp(['Variable ', u_varname, ' already exists. Updating its value']);
end

if (~got.Vvar)
  ic = ic + 1;
  S.Dimensions = I.Dimensions;
  S.Variables(ic) = roms_metadata(v_varname, spherical);
  append_vars = true;
else
  disp(['Variable ', v_varname, ' already exists. Updating its value']);
end

if (append_vars)
  check_metadata(S);
  nc_append(ncfile, S);
end

%--------------------------------------------------------------------------
%  Write out A-grid vector component into NetCDF file.
%--------------------------------------------------------------------------

%  Compute vector components at RHO-points and rotate to Eastward and
%  Nortward directions

Nrec = double(size(nc_read(ncfile, 'ocean_time')));

for rec = 1:Nrec

  switch u_varname
    case 'u_eastward'
      Uinp = nc_read(ncfile, 'u', rec);
    case 'ubar_eastward'
      Uinp = nc_read(ncfile, 'ubar', rec);
  end

  switch v_varname
    case 'v_northward'
      Vinp = nc_read(ncfile, 'v', rec);
    case 'vbar_northward'
      Vinp = nc_read(ncfile, 'vbar', rec);
  end

  if (length(size(Uinp)) == 3)
    is3d = true;
    N = size(Uinp,3);
    angle = repmat(G.angle, [1,1,N]);
  else
    is3d = false;
    angle = G.angle;
  end
  SinAngle = sin(angle);
  CosAngle = cos(angle);
  
  Urho = NaN(size(angle));
  Vrho = NaN(size(angle));

  [Lp,Mp] = size(G.angle);
  L = Lp-1;  Lm = L-1;
  M = Mp-1;  Mm = M-1;

  % Average to RHO-points, A-grid. Apply gradient lateral boundary.
  
  if (is3d)
    Urho(2:L, 1:Mp, :) = 0.5 .* (Uinp(1:Lm, 1:Mp, :)+Uinp(2:L,  1:Mp, :));
    Urho(1,   :,    :) = Urho(2, :, :);
    Urho(Lp,  :,    :) = Urho(L, :, :);

    Vrho(1:Lp, 2:M, :) = 0.5 .* (Vinp(1:Lp, 1:Mm, :)+Vinp(1:Lp, 2:M,  :));
    Vrho(:,    1,   :) = Vrho(:, 2, :);   
    Vrho(:,    Mp,  :) = Vrho(:, M, :);   
  else
    Urho(2:L,  1:Mp) = 0.5 .* (Uinp(1:Lm, 1:Mp)+Uinp(2:L,  1:Mp));
    Urho(1,   :) = Urho(2, :);
    Urho(Lp,  :) = Urho(L, :);

    Vrho(1:Lp, 2:M ) = 0.5 .* (Vinp(1:Lp, 1:Mm)+Vinp(1:Lp, 2:M ));
    Vrho(:,    1 ) = Vrho(:, 2);   
    Vrho(:,    Mp) = Vrho(:, M);   
  end
  
  % Rotate geographical Eastward and Northward components.

  Urot = Urho .* CosAngle - Vrho .* SinAngle;
  Vrot = Vrho .* CosAngle + Urho .* SinAngle;
  
  % Write out rotated A-grid components.

  nc_write(ncfile, u_varname, Urot, rec);
  nc_write(ncfile, v_varname, Vrot, rec);

  F.time = nc_read(ncfile, 'ocean_time', rec);
  F(rec).Urho = Urho;
  F(rec).Vrho = Vrho;
  F(rec).Urot = Urot;
  F(rec).Vrot = Vrot;
end

return
