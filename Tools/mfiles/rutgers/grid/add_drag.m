function add_drag(ncfile, D)

%
% ADD_DRAG:  Adds bottom drag variables to a ROMS Grid NetCDF file
%
% add_drag(ncfile, D)
%
% This function adds time invariant, spatially-varying bottom drag 
% variable(s) (ZoBot, rdrag, or rdrag2) to an existing ROMS Grid NetCDF
% file. Any of these variables are used when UV_DRAG_GRID is activated
% in ROMS. The bottom roughness length (ZoBot), linear bottom drag
% coefficients (rdrag), or quadratic bottom drag coefficients (rdrag2)
% is needed if the option UV_LOGDRAG, UV_LDRAG, or UV_QDRAG is activated,
% respectively.
%
% This function will add all variables available in input structure, D.
%
% On Input:
%
%    ncfile        GRID NetCDF file name (string)
%
%    D             Bottom drag variables at RHO-points (struct):
%
%                    D.ZoBot     bottom roughness length (m)
%                    D.rdrag     linear bottom drag coefficients (m/s)
%                    D.rdrag2    quadratic bottom drag coefficients
%

% svn $Id: add_drag.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire grid NetCDF file.
%--------------------------------------------------------------------------

I = nc_inq(ncfile);

got.ZoBot  = any(strcmp({I.Variables.Name}, 'ZoBot'));
got.rdrag  = any(strcmp({I.Variables.Name}, 'rdrag'));
got.rdrag2 = any(strcmp({I.Variables.Name}, 'rdrag2'));

if (any(strcmp({I.Variables.Name}, 'spherical'))),
  spherical = nc_read(ncfile, 'spherical');
  if (ischar(spherical)),
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

if isfield(D, 'ZoBot'),
  if (~got.ZoBot),
    ic = ic + 1;
    S.Dimensions = I.Dimensions;
    S.Variables(ic) = roms_metadata('ZoBot', spherical);
    append_vars = true;
  else
    disp(['Variable "ZoBot" already exists. Updating its value']);
  end
end

if isfield(D, 'rdrag'),
  if (~got.rdrag),
    ic = ic + 1;
    S.Dimensions = I.Dimensions;
    S.Variables(ic) = roms_metadata('rdrag', spherical);
    append_vars = true;
  else
    disp(['Variable "rdrag" already exists. Updating its value']);
  end
end

if isfield(D, 'rdrag2'),
  if (~got.rdrag),
    ic = ic + 1;
    S.Dimensions = I.Dimensions;
    S.Variables(ic) = roms_metadata('rdrag2', spherical);
    append_vars = true;
  else
    disp(['Variable "rdrag2" already exists. Updating its value']);
  end
end

if (append_vars),
  check_metadata(S);
  nc_append(ncfile, S);
end

%--------------------------------------------------------------------------
%  Write out sponge variables into GRID NetCDF file.
%--------------------------------------------------------------------------

Lr = I.Dimensions(strcmp({I.Dimensions.Name},'xi_rho' )).Length;
Mr = I.Dimensions(strcmp({I.Dimensions.Name},'eta_rho')).Length;

if isfield(D, 'ZoBot'),
 [Im,Jm] = size(D.ZoBot);
 if (Im == Lr && Jm == Mr),
    nc_write(ncfile, 'ZoBot', D.ZoBot);
 else
   error([' ADD_DRAG: size(D.ZoBot) is different to Lr = ',             ...
         num2str(Lr), blanks(3), 'Mr = ', num2str(Mr)]);	 
 end
end

if isfield(D, 'rdrag'),
 [Im,Jm] = size(D.rdrag);
 if (Im == Lr && Jm == Mr),
    nc_write(ncfile, 'rdrag', D.rdrag);
 else
   error([' ADD_DRAG: size(D.rdrag) is different to Lr = ',             ...
         num2str(Lr), blanks(3), 'Mr = ', num2str(Mr)]);	 
 end
end

if isfield(D, 'rdrag2'),
 [Im,Jm] = size(D.rdrag2);
 if (Im == Lr && Jm == Mr),
    nc_write(ncfile, 'rdrag2', D.rdrag2);
 else
   error([' ADD_DRAG: size(D.rdrag2) is different to Lr = ',            ...
         num2str(Lr), blanks(3), 'Mr = ', num2str(Mr)]);	 
 end
end

return
