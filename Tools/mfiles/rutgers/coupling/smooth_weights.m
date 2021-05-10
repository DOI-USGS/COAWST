function [Dweights, Oweights] = smooth_weights (Wname, varargin)

%
% SMOOTH_WEIGHTS: Sets smooth weight factors for merging DATA and ROMS
%                 components
%
% [Dweights, Oweights] = smooth_weights(Wname, Nconv, Lplot)
%
% This function computes smooth melding weights to combine fields from
% DATA and OCEAN components. The merging factors change in an area next
% to the OCEAN grid open boundaries for a smooth transition.
%
% WARNING: The output "Dweights" and "Oweights" are not written into the
%          weights NetCDF file. It is done elsewhere.
%
% On Input:
%
%    Wname       Weights NetCDF filename created elsewhere (string)
%            or, an existing weights structure computed elsewhere (struc)
%
%    Nconv       Number of smoothing convolutions iterations (Optional)
%                  (default = 50; higher values for wider melding zone)
%
%    Lplot       Switch to plot melding weights (Optional)
%                  (default = true)
%
% On Output:
%
%    Dweights    DATA  component smooth melding weights (2D array)
%
%    Oweights    OCEAN component smooth melding weights (2D array)
%  

% svn $Id: smooth_weights.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           John L. Wilkin        %
%=========================================================================%

% Initialize.

switch numel(varargin)
  case 0
    Nconv = 50;
    Lplot = true;
  case 1
    Nconv = varargin{1};
    Lplot = true;
  case 2
    Nconv = varargin{1};
    Lplot = varargin{2};
end

if (~isstruct(Wname)),
  List = {'lon', 'lat', 'mask', 'data_weight_rigid', 'ocean_weight_rigid'};
  for var = List
    field = char(var);
    W.(field) = ncread(Wname, field);
  end
else
  W = Wname;
end

%--------------------------------------------------------------------------
% Compute smooth melding weights.
%--------------------------------------------------------------------------

% Sum atmosphere mask (0: land, 1:sea) and DATA rigid weigths. 

M = W.mask + W.data_weight_rigid;         % 0 = land; 1 = ROMS; 2 = data

S = M;  

for k=1:Nconv
  
  S = conv2(S, ones([3 3])/9, 'same');
 
  % Keep points with 1<S<2 : they are on the ROMS perimeter/data boundary
  % All others restore before next pass of filter
  
  % restore ROMS coastal mask

  tmp = find(S>0&S<1);
  S(tmp) = M(tmp);
  
  % restore land and data mask

  S(M==2) = 2;
  S(M==0) = 0;
  
  % restore perimeter

  S(:,[1 end]) = M(:,[1 end]);
  S([1 end],:) = M([1 end],:);
  
end

% Shift so that 0<S<1

S = max(0,S-1);

% Update weights so that melding zone inside ROMS perimeter adjacent to
% data gets the values of S>0 and S<1. Values sum to 1. 

% Smoothed data weights

data_weight_smooth = W.data_weight_rigid;
melded = (W.data_weight_rigid==0);
data_weight_smooth(melded) = S(melded);

% Smoothed ocean weights

ocean_weight_smooth = W.ocean_weight_rigid;
melded = (W.ocean_weight_rigid==1);
ocean_weight_smooth(melded) = 1-S(melded);

%--------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------

if (Lplot)
  figure;
  
  hax(1) = subplot(2,2,1);
  pcolorjw(W.lon, W.lat, W.ocean_weight_rigid); colorbar
  title('ocean\_weight\_rigid')

  hax(2) = subplot(2,2,2);
  pcolorjw(W.lon, W.lat, ocean_weight_smooth); colorbar
  title(['ocean\_weight\_smooth, N = ', num2str(Nconv)]);
  
  hax(3) = subplot(2,2,3);
  pcolorjw(W.lon, W.lat, W.data_weight_rigid); colorbar
  title('data\_weight\_rigid')

  hax(4) = subplot(2,2,4);
  pcolorjw(W.lon, W.lat, data_weight_smooth); colorbar
  title(['data\_weight\_smooth, N = ', num2str(Nconv)]);

  linkaxes(hax)
end

%--------------------------------------------------------------------------
% Set output arguments.
%--------------------------------------------------------------------------

Dweights = data_weight_smooth;
Oweights = ocean_weight_smooth;

return
