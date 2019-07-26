function [vint,varargout] = roms_zint_mw(data,grd,varargin)
% $Id$
% ROMS_ZINT:  vertical integral at rho points from depth z1 to z2
%
% USAGE:  [vint,thickness,dz] = roms_zint212(data,grd);
%
%         In this case, the vertical integral is over the entire
%         water column.
%
% USAGE:  [vint,thickness,dz] = roms_zint212(data,grd,z_r_indices);
% 
%         In this case, the vertical integral is between two sigma
%         levels.
%
% USAGE:  [vint,thickness,dz] = roms_zint212(data,grd,z2,z1);
%
%         In this case, the vertical integral is between specified
%         depths.  They can either be two scalar values or two general
%         arrays.  
%
% PARAMETERS:
% Input:
%     data: 3D variable from ROMS file
%     grd:  as returned from roms_get_grid.m
%     z_r_indices:
%         Optional.  If provided, then this is a tuple with the indices
%         into the grd.z_r array.  We wish to integrate between the 
%         two levels specified by these indices.
%     z2, z1:
%         Optional.  If not supplied, then the integral is over the
%         whole water column.
%        
%         If either z1 or z2 are arrays, then the output argument dz
%         cannot be asked for.
%
%         If z1 or z2 are 2-d arrays is it assumed they described a 
%         matrix of depths (negative) defined on the rho points grid, 
%         such as the depth of some property (mixed layer, an isopycnal 
%         surface) determined elsewhere.  z1 could be zeta and z2 could 
%         be -h, in which we just get the integral over the full water 
%         column.
%
%         If z1 == 0 or NaN it is assumed the integral is from the 
%         surface to z2.
%
%         If z2 == NaN or [] it is assumed the integral is from z1 to 
%         -h.
%
%        The order of integration is from z2 to z1, so please make z2
%        the third argument, and z1 the fourth.
% Output:
%     vint:
%         Vertical integral.
%     thickness:
%         Optional.  2D array of water column heights between the given 
%         limits of integration.  This allows a vertical average to be 
%         calculated if desired.
%     dz:
%         Thicknesses of each layer in the context of z_w.  We cannot
%         compute this if z1 and z2 are given as arrays.
%
%
% Method is to compute a matrix of layer thicknesses, then multiply the
% data element-by-element, then sum over the vertical dimension
%
% To integrate over depths not the full water column, the layer thickness
% matrix is modified accordingly by effectively implementing
%
% int_z2^z1 dz = int_-h^0 dz - int_z1^0 dz - int_-h^z2 dz

%
% Set up defaults.
INTEGRATE_OVER_ENTIRE_WATER_COLUMN = 0;  % do not integrate over entire water column.
INTEGRATE_BETWEEN_LEVELS = 0;            % Do not assume that the limits of integration
                                         % line up nicely with z_r levels.

%
% check the inputs.  Must have at least two.
if (nargin < 2 )
  msg = sprintf ( '%s:  not enough input arguments.\n', mfilename );
  error ( msg );
end

%
% If this case, integrate over the entire water column.  Set the z_r indices
% to cover the entire water column.  We can also compute dz in this case.
if nargin == 2
  INTEGRATE_OVER_ENTIRE_WATER_COLUMN = 1;
end

%
% If this case, then we were given the indices into z_r.
% We can still compute dz since we are given indices into the vertical
% laters of z_r.
if ( nargin == 3 )
  z_r_indices = varargin{1};
  if ( length(z_r_indices) ~= 2 )
    msg = sprintf ( '%s:  length of third argument must be two.\n', mfilename );
    error ( msg );
  end
  z_r_index_2 = z_r_indices(1);
  z_r_index_1 = z_r_indices(2);
  
  %
  % check that the indices make sense
  if z_r_index_2 > z_r_index_1
    msg = sprintf ( '%s:  index 2 (%d) must be less than index 1 (%d).\n', mfilename, z_r_index_2, z_r_index_1 );
    error ( msg );
  end
  INTEGRATE_BETWEEN_LEVELS = 1;
end


%
% Check the position of the data on the Arakawa c-grid
pos = roms_cgridpos(data,grd);
switch pos
  case 'u'
    z_w = 0.5*(grd.z_w(:,1:end-1,:)+grd.z_w(:,2:end,:));
    z_r = 0.5*(grd.z_r(:,1:end-1,:)+grd.z_r(:,2:end,:));
    h = 0.5*(grd.h(1:end-1,:)+grd.h(2:end,:));
  case 'v'
    z_w = 0.5*(grd.z_w(:,:,1:end-1)+grd.z_w(:,:,2:end));
    z_r = 0.5*(grd.z_r(:,:,1:end-1)+grd.z_r(:,:,2:end));
    h = 0.5*(grd.h(:,1:end-1)+grd.h(:,2:end));
  otherwise
    z_w = grd.z_w;
    z_r = grd.z_r;
    h = grd.h;
end

%
% In this case, we were given two 2D arrays, between which we wish to
% integrate.
if ( nargin == 4 )
  
  z2 = varargin{1};
  z1 = varargin{2};
  
  %
  % check that both are 2D
  if ( (ndims(z1) ~= 2) | (ndims(z2) ~= 2 ) )
    msg = sprintf ( '%s:  3rd and 4th arguments cannot be 3D or higher.\n', mfilename );
    error ( msg );
  end
  
  %
  % Get the sizes of the z levels.
  % [tz,rz,cz] = size(grd.z_r);
  [tz,rz,cz] = size(z_r);
  [rz1,cz1] = size(z1);
  [rz2,cz2] = size(z2);
  
  %
  % Check if we were given scalar values for z1 and z2
  if ( ( rz1 == 1 ) & (cz1 == 1) & (rz2 == 1) & (cz2 == 1) )
    z1 = z1*ones(rz,cz);
    z2 = z2*ones(rz,cz);
    [rz1,cz1] = size(z1);
    [rz2,cz2] = size(z2);
  end
  
  %
  % Check that X and Y sizes are the same
  if ( ( rz1 ~= rz ) | ( cz1 ~= cz ) )
    format = '%s:  size of z1 [%dx%d] does not match YX dimensions of z_r [%dx%dx%d]\n';
    msg = sprintf ( format, mfilename, rz1, cz1, tz, rz, cz );
    error ( msg );
  end
  if ( ( rz2 ~= rz ) | ( cz2 ~= cz ) )
    format = '%s:  size of z2 [%dx%d] does not match YX dimensions of z_r [%dx%dx%d]\n';
    msg = sprintf ( format, mfilename, rz2, cz2, tz, rz, cz );
    error ( msg );
  end
  
  
  
end

%
% No more than four input arguments.
if ( nargin > 4 )
  msg = sprintf ( '%s:  too many input arguments.\n', mfilename );
  error ( msg );
end


%
% Must have at least one output argument
if (nargout < 1)
  msg = sprintf ( '%s:  must have at least 1 output argument.\n', mfilename );
  error ( msg );
end
if (nargout > 3)
  msg = sprintf ( '%s:  too many output arguments.\n', mfilename );
  error ( msg );
end



% int_-h^0 dz first:
% layer thicknesses
%grd.z_w(1,:,:) = -65 * ones(130,162);
%grd.z_w(2,:,:) = -64 * ones(130,162);
%grd.z_w(3,:,:) = -63 * ones(130,162);
%grd.z_w(4,:,:) = -61 * ones(130,162);
%grd.z_w(5,:,:) = -59 * ones(130,162);
%grd.z_w(6,:,:) = -56 * ones(130,162);
%grd.z_w(7,:,:) = -52 * ones(130,162);
%grd.z_w(8,:,:) = -50 * ones(130,162);
%grd.z_w(9,:,:) = -48 * ones(130,162);
%grd.z_w(10,:,:) = -45 * ones(130,162);
%grd.z_w(11,:,:) = -42 * ones(130,162);
%grd.z_w(12,:,:) = -32 * ones(130,162);
%grd.z_w(13,:,:) = -22 * ones(130,162);
%grd.z_w(14,:,:) = -12 * ones(130,162);
%grd.z_w(15,:,:) = -11 * ones(130,162);
%grd.z_w(16,:,:) = -10 * ones(130,162);
%grd.z_w(17,:,:) = -9 * ones(130,162);
%grd.z_w(18,:,:) = -8 * ones(130,162);
%grd.z_w(19,:,:) = -7 * ones(130,162);
%grd.z_w(20,:,:) = -5 * ones(130,162);
%grd.z_w(21,:,:) = -1 * ones(130,162);


% dz = diff(grd.z_w,1,1); 
dz = diff(z_w,1,1); 

%	data(1,:,:) = 3 * ones(130, 162);
%	data(2,:,:) = 4 * ones(130, 162);
%	data(3,:,:) = 5 * ones(130, 162);
%	data(4,:,:) = 6 * ones(130, 162);
%	data(5,:,:) = 7 * ones(130, 162);
%	data(6,:,:) = 8 * ones(130, 162);
%	data(7:20,:,:) = ones(14,130,162);



%
% Compute the vertical integral int_z1^z2 over the entire water column.
vint = data.*dz;
vint = squeeze(sum(vint,1));




%
% Now adjust the value of the vertical integral depending upon the 
% limits of integration.
if INTEGRATE_OVER_ENTIRE_WATER_COLUMN
  
  thickness = squeeze(sum(dz,1));
  
  %
  % We're basically done.  Just pack up the output arguments and
  % be done with it.
  switch nargout
    case 1
      % 
      % do nothing
      ;
    case 2
      varargout{1} = thickness;
    case 3
      
      % thickness is sum of dz
      varargout{1} = thickness;
      varargout{2} = dz;
      
  end
  
  return
end


if INTEGRATE_BETWEEN_LEVELS
  
  
  [nr,r,c] = size(data);
  %
  % subtract int_z1^0 and int_-h^z2 and we're done
  %
  % If the z_r_index_1 is N, then don't bother.  
  % It just produces a 0-row matrix.
  %
  % The notation is a little confusing, I admit.
  if ( z_r_index_1 == grd.N )
    int_z1_0 = zeros(r,c);
  else
    dz_z1_0 = dz(z_r_index_1+1:end,:,:);
    int_z1_0 = data(z_r_index_1+1:end,:,:) .* dz_z1_0;
    int_z1_0 = squeeze(sum(int_z1_0,1));
  end
  
  % If the z_r_index_2 is 1, then don't bother.  
  % It just produces a 0-row matrix.
  if ( z_r_index_2 == 1 )
    int_h_z2 = zeros(r,c);
  else
    dz_h_z2 = dz(1:z_r_index_2-1,:,:);
    int_h_z2 = data(1:z_r_index_2-1,:,:) .* dz_h_z2;
    int_h_z2 = squeeze(sum(int_h_z2,1));
  end
  
  
  vint = vint - int_z1_0 - int_h_z2;
  
  %
  % dz is restricted to a certain subset of levels.
  dz = dz(z_r_index_2:z_r_index_1,:,:);
  thickness = squeeze(sum(dz,1));
  
  
  %
  % We're basically done.  Just pack up the output arguments and
  % be done with it.
  switch nargout
    case 1
      % 
      % do nothing
      ;
    case 2
      varargout{1} = thickness;
    case 3
      
      % thickness is sum of dz
      varargout{1} = thickness;
      varargout{2} = dz;
      
  end
  
  return
end


%
% We find the z_w values that bound the limits of integration.  Linearly
% interpolate to find how much of the bounding dz value lies between the
% bounding z_w values and adjust the dz value by that much.

% num_depths = size(grd.z_r,1);
num_depths = size(z_r,1);

real_bottom_layer = z2;
real_top_layer = z1;

%
% find the indicies that bracket z2
alpha = [1 0]';
for r = 1:rz
  for c = 1:cz
    
    %
    % find the indicies that bracket z2
    ind = find(z_w(:,r,c) <= z2(r,c) );
    if isempty(ind)
      %
      % do nothing, the entire water column is bounded below
      % by z2.  
      % But, thickness makes less sense here
      % since there's a gap between z2 and the bottom level.
      % So adjust thickness here from the default of z1-z2;
      real_bottom_layer(r,c) = -1*h(r,c);
      ; 
    elseif ( length(ind) > grd.N )
      %
      % z2 lies above the water column.  No contribution at all
      % is possible here.
      dz(:,r,c) = 0;
    else
      
      %
      % linearly adjust the dz value that's bracketed
      % Values below that dz have no contribution.
      bottom = ind(end);
      top = ind(end) + 1;
      br_inds = [bottom top];
      z_br = z_w(br_inds,r,c);
      alpha_i = interp1 ( z_br, alpha, z2(r,c) );
      dz(bottom,r,c) = dz(bottom,r,c) * alpha_i;
      dz(1:bottom-1,r,c) = 0;
    end
    
    
    %
    % find the indicies that bracket z1, which lies above z2
    ind = find(z_w(:,r,c) > z1(r,c) );
    if isempty(ind)
      %
      % do nothing, the entire water column is bounded above
      % by z1.  That's ok.
      % But, thickness makes less sense here
      % since there's a gap between z2 and the bottom level.
      % So adjust thickness here from the default of z1-z2;
      real_top_layer(r,c) = 0;
      ; 
    elseif ( length(ind) > grd.N )
      %
      % z1 lies below the water column.  No contribution at all
      % is possible here.
      dz(:,r,c) = 0;
    else
      
      %
      % linearly adjust the dz value that's bracketed
      % Values above that dz have no contribution.
      bottom = ind(1)-1;
      top = ind(1);
      br_inds = [bottom top];
      z_br = z_w(br_inds,r,c);
      alpha_i = interp1 ( z_br, alpha, z1(r,c) );
      dz(bottom,r,c) = dz(bottom,r,c) * (1 - alpha_i);
      dz(bottom+1:end,r,c) = 0;
    end
    
  end
end


%
% Just redo the integral now.
vint = data.*dz;
vint = squeeze(sum(vint,1));
thickness = real_top_layer - real_bottom_layer;


if nargout == 2
  varargout{1} = thickness;
elseif nargout == 3
  varargout{2} = dz;
end

return;

