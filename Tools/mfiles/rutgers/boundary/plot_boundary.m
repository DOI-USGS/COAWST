function plot_boundary(Gname, Bname, Vname, Tindex)

%
% PLOT_BOUNDARY:  Plots requested ROMS variable from boundary NetCDF file
%
% B=plot_boundary(Gname, Bname, Vname, Tindex)
%
% This function plots requested ROMS variable from input lateral
% boundary NetCDF file. This function is very useful for browsing
% 3D lateral boundary conditions for ROMS in a vertical slab. The
% location of the minumum is marked with a filled magenta circle
% whereas the  maximum is marked with a filled magenta square.
%
% On Input:
%
%    Gname         ROMS Grid NetCDF file/URL name (string)
%              or, an existing ROMS grid structure (struct array)
%
%    Bname         ROMS boundary NetCDF file/URL name (string)
%
%    Vname         ROMS state variable name to process (string)
%
%    Tindex        Time record index used to process (scalar)
%


% svn $Id: plot_boundary.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

Tname = [];
Tsize = 0;

% Set ROMS Grid structure.
  
if (~isstruct(Gname))
  G = get_roms_grid(Gname,Bname);
else
  G = Gname;
end

% Get boundary NetCDF file information.

I = nc_inq(Bname);

% Determine available boundary variables for requested state variable.

Vnames = {I.Variables.Name};
Vindex = find(strncmp(Vnames, strcat(Vname,'_'), length(Vname)+1));

if (~isempty(Vindex))
  Bvars = Vnames(Vindex);
  Nvars = length(Bvars);
else
  error(['PLOT_BOUNDARY: unable to find boundary variables for: ',Vname]);
end

%--------------------------------------------------------------------------
% Read and plot boundary requested lateral boundary variable(s).
%--------------------------------------------------------------------------

for n = 1:Nvars
   
  VarName  = char(Bvars(n));
  index    = Vindex(n);
  nvdims   = length(I.Variables(index).Dimensions);
  is3d     = false;

  western  = contains(Bvars(n), '_west' );
  eastern  = contains(Bvars(n), '_east' );
  southern = contains(Bvars(n), '_south');
  northern = contains(Bvars(n), '_south');
  
  if (nvdims > 0)
    for m = 1:nvdims
      dname = char(I.Variables(index).Dimensions(m).Name);
      switch dname
        case {'s_rho'}
          is3d = true;
          Km   = I.Variables(index).Dimensions(m).Length;
        case {'xi_rho'}
          if (G.spherical)
            if (southern)
              XorY = G.lon_rho_south;
            elseif (northern)
              XorY = G.lon_rho_north;
            end
          else
            if (southern)
              XorY = G.x_rho_south;
            elseif (northern)
              XorY = G.x_rho_north;
            end
          end
          if (southern)
            mask = squeeze(G.mask_rho(:,1));
            Z    = squeeze(G.z_r(:,1,:));
          elseif (northern)
            mask = squeeze(G.mask_rho(:,end));
            Z    = squeeze(G.z_r(:,end,:));
          end
        case {'eta_rho'}
          if (G.spherical)
            if (western)
              XorY = G.lat_rho_west';
            else
              XorY = G.lat_rho_east';
            end
          else
            if (western)
              XorY = G.y_rho_west';
            else
              XorY = G.y_rho_east';
            end
          end
          if (western)
            mask = squeeze(G.mask_rho(1,:));
            Z    = squeeze(G.z_r(1,:,:));
          elseif (eastern)
            mask = squeeze(G.mask_rho(end,:));
            Z    = squeeze(G.z_r(end,:,:));
          end
        case {'xi_u'}
          if (G.spherical)
            if (southern)
              XorY = G.lon_u_south;
            elseif (northern)
              XorY = G.lon_u_north;
            end
          else
            if (southern)
              XorY = G.x_u_south;
            elseif (northern)
              XorY = G.x_u_north;
            end
          end
          if (southern)
            mask = squeeze(G.mask_u(:,1));
            Z    = squeeze(G.z_u(:,1,:));
          elseif (northern)
            mask = squeeze(G.mask_u(:,end));
            Z    = squeeze(G.z_u(:,end,:));
          end
        case {'eta_u'}
          if (G.spherical)
            if (western)
              XorY = G.lat_u_west';
            elseif (eastern)
              XorY = G.lat_u_east';
            end
          else
            if (western)
              XorY = G.y_u_west';
            elseif (eastern)
              XorY = G.y_u_east';
            end
          end
          if (western)
            mask = squeeze(G.mask_u(1,:));
            Z    = squeeze(G.z_u(1,:,:));
          elseif (eastern)
            mask = squeeze(G.mask_u(end,:));
            Z    = squeeze(G.z_u(end,:,:));
          end
        case {'xi_v'}
          if (G.spherical)
            if (southern)
              XorY = G.lon_v_south;
            elseif (northern)
              XorY = G.lon_v_north;
            end
          else
            if (southern)
              XorY = G.x_v_south;
            elseif (northern)
              XorY = G.x_v_north;
            end
          end
          if (southern)
            mask = squeeze(G.mask_v(:,1));
            Z    = squeeze(G.z_v(:,1,:));
          elseif (northern)
            mask = squeeze(G.mask_v(:,end));
            Z    = squeeze(G.z_v(:,end,:));
          end
        case {'eta_v'}
          if (G.spherical)
            if (western)
              XorY = G.lat_v_west;
            elseif (eastern)
              XorY = G.lat_v_east;
            end
          else
            if (western)
              XorY = G.y_v_west';
            elseif (eastern)
              XorY = G.y_v_east';
            end
          end
          if (western)
            mask = squeeze(G.mask_v(1,:));
            Z    = squeeze(G.z_v(1,:,:));
          elseif (eastern)
            mask = squeeze(G.mask_v(end,:));
            Z    = squeeze(G.z_v(end,:,:));
          end
        case 'bry_time'
          Tsize = I.Variables(index).Dimensions(m).Length;
          Tname = dname;
          Tind  = strcmp(Vnames,Tname);
          Aind  = strcmp({I.Variables(Tind).Attributes.Name},'units');
          Aval  = I.Variables(Tind).Attributes(Aind).Value;
          if (contains(Aval, 'since'))
            epoch = datenum(Aval(strfind(Aval,'since')+6:end));
          else
            epoch = 0;
          end
      end
    end
  end

% Read in lateral boundary variable.

  if (Tindex > Tsize)
    Tindex = Tsize;                     % process last time record available
  end 

  ReplaceValue = NaN;
  PreserveType = false;

  B    = nc_read(Bname,VarName,Tindex,ReplaceValue,PreserveType);
  time = nc_read(Bname,Tname,Tindex);

  if (epoch > 0)
    my_date = datestr(epoch+time/86400);
  else
    my_date = num2str(time/86400);
  end
    
  Bmin = min(B(:));  [Imin,Jmin] = find(B == min(B(:))); 
  Bmax = max(B(:));  [Imax,Jmax] = find(B == max(B(:)));

  if (length(Imin) > 1 && length(Jmin > 1))
    Imin = Imin(1);
    Jmin = Jmin(1);
  end

  if (length(Imax) > 1 && length(Jmax > 1))
    Imax = Imax(1);
    Jmax = Jmax(1);
  end

  if (is3d)
    mask = repmat(mask(:), [1,Km]);   
    XorY = repmat(XorY(:), [1,Km]);   
  end

  ind = find(mask < 0.5);
  if (~isempty(ind))
    B(ind) = NaN;
  end
  
% Plot lateral boundary data.

  if (is3d)
    figure;
    pcolor(XorY,Z,B); shading interp; colorbar

    ht = title([VarName, ':', blanks(4),                                ...
                'Record = ', num2str(Tindex), ',', blanks(4),           ...
                my_date],                                               ...
                'FontSize', 14, 'FontWeight', 'bold' );  
    set(ht, 'interpreter', 'none');                    % no TeX interpreter

    hx = xlabel(['Min = ', num2str(Bmin), blanks(4),                    ...
                 '(',num2str(Imin),', ',num2str(Jmin), '),', blanks(8), ...
                 'Max = ', num2str(Bmax), blanks(4),                    ...
                 '(', num2str(Imax), ', ', num2str(Jmax), ')'],         ...
                 'FontSize', 14, 'FontWeight', 'bold' );
    
    hold on;
    hm = plot(XorY(Imin,Jmin), Z(Imin,Jmin), 'o',                       ...
              XorY(Imax,Jmax), Z(Imax,Jmax), 's');
    set(hm,'MarkerEdgeColor','none','MarkerFaceColor','m','MarkerSize',12);

    hb = plot(XorY(Imin,Jmin), Z(Imin,Jmin), 'ko',                      ...
              XorY(Imax,Jmax), Z(Imax,Jmax), 'ks');
    set(hb,'MarkerSize',12);  
    hold off;
  
  else

    figure;
    plot(XorY, B, 'r-')

    ht = title([VarName, ':', blanks(4),                                ...
                'Record = ', num2str(Tindex), ',', blanks(4),           ...
                datestr(epoch+time/86400)],                             ...
                'FontSize', 14, 'FontWeight', 'bold' );  
    set(ht, 'interpreter', 'none');                    % no TeX interpreter
  
  end

end

return
