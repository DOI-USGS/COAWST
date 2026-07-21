function [nfast,weight] = set_weights(ndtfast, varargin)

% SET_WEIGHTS:  ROMS weights for barotropic time averaged fields
%
% [nfast,weight] = set_weights(ndtfast, Lplot)
%
%  This function computes the barotropic kernel time averaging weight
%  using the power-law shape filters, which are computed as:
%
%     F(xi)=xi^Falpha*(1-xi^Fbeta)-Fgamma*xi
%
%  where xi=scale*i/ndtfast; and scale, Falpha, Fbeta, Fgamma, and
%  normalization are chosen to yield the correct zeroth-order
%  (normalization), first-order (consistency), and second-order moments,
%  resulting in overall second-order temporal accuracy for time-averaged
%  barotropic motions resolved by baroclinic time step.
%
% On Input:
%  
%    ndtfast     Number of barotropic substeps for each baroclinic step
%                  (integer)
%    Lplot       Switch to plot the weights (optional)
%                  default = false
%
% On Output:
%
%    nfast       Number of baroclinic substeps needed for 2D kernel
%                  sveraging (integer)
%    weights     Normalized averaging weights (real array)
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%  
  
%  Initialize.

switch numel(varargin)
  case 0
    Lplot=false;
  case 1
    Lplot=varargin{1};
end  
  
Falpha=2.0;
Fbeta =4.0;
Fgamma=0.284;

scale=(Falpha+1.0)*(Falpha+Fbeta+1.0)/                                  ...
      ((Falpha+2.0)*(Falpha+Fbeta+2.0)*double(ndtfast));

%  Find center of gravity of the primary weighting shape function and
%  iteratively adjust "scale" to place the  centroid exactly at
%  "ndtfast".

%gamma=Fgamma*max(0.0, 1.0-10.0/double(ndtfast));
 gamma=Fgamma*(1.0-2.8/double(ndtfast));

weight=zeros([2 256]);
 
for iter=1:16
  nfast=0;
  for i=1:2*ndtfast
    cff=scale*double(i);
    weight(1,i)=cff^Falpha-cff^(Falpha+Fbeta)-gamma*cff;
    if (weight(1,i) > 0.0)
      nfast=i;
    end
    if ((nfast > 0) && (weight(1,i) < 0.0))
      weight(1,i)=0.0;
    end
  end
  wsum=0.0;
  shift=0.0;
  for i=1:nfast
    wsum=wsum+weight(1,i);
    shift=shift+weight(1,i)*double(i);
  end  
  scale=scale*shift/(wsum*double(ndtfast));
end

%-----------------------------------------------------------------------
%  Post-processing of primary weights.
%-----------------------------------------------------------------------
%
%  Although it is assumed that the initial settings of the primary
%  weights has its center of gravity "reasonably close" to NDTFAST,
%  it may be not so according to the discrete rules of integration.
%  The following procedure is designed to put the center of gravity
%  exactly to NDTFAST by computing mismatch (NDTFAST-shift) and
%  applying basically an upstream advection of weights to eliminate
%  the mismatch iteratively. Once this procedure is complete primary
%  weights are normalized.
%
%  Find center of gravity of the primary weights and subsequently
%  calculate the mismatch to be compensated.

for iter=1:ndtfast
  wsum=0.0;
  shift=0.0;
  for i=1:nfast
    wsum=wsum+weight(1,i);
    shift=shift+double(i)*weight(1,i);
  end
  shift=shift/wsum;
  cff=double(ndtfast)-shift;

%  Apply advection step using either whole, or fractional shifts.
%  Notice that none of the four loops here is reversible.

  if (cff > 1.0)
    nfast=nfast+1;
    for i=nfast:-1:2
      weight(1,i)=weight(1,i-1);
    end
    weight(1,1)=0.0;
  elseif (cff > 0.0)
    wsum=1.0-cff;
    for i=nfast:-1:2
      weight(1,i)=wsum*weight(1,i)+cff*weight(1,i-1);
    end
    weight(1,1)=wsum*weight(1,1);
  elseif (cff < -1.0)
    nfast=nfast-1;
    for i=1:nfast
      weight(1,i)=weight(1,i+1);
    end
    weight(1,nfast+1)=0.0;
  elseif (cff < 0.0)
    wsum=1.0+cff;
    for i=1:nfast-1
      weight(1,i)=wsum*weight(1,i)-cff*weight(1,i+1);
    end
    weight(1,nfast)=wsum*weight(1,nfast);
  end
end

disp(['  ndtfast = ', num2str(ndtfast), blanks(3),                      ...
      'nfast = ', num2str(nfast)]);

%  Set SECONDARY weights assuming that backward Euler time step is used
%  for free surface.  Notice that array weight(2,i) is assumed to
%  have all-zero status at entry in this segment of code.

for j=1:nfast
  cff=weight(1,j);
  for i=1:j
    weight(2,i)=weight(2,i)+cff;
  end
end

%  Normalize both set of weights.

wsum=0.0;
cff=0.0;
for i=1:nfast
  wsum=wsum+weight(1,i);
  cff=cff+weight(2,i);
end
wsum=1.0/wsum;
cff=1.0/cff;
for i=1:nfast
  weight(1,i)=wsum*weight(1,i);
  weight(2,i)=cff*weight(2,i);
end

% Plot weights.

if (Lplot)
  figure;

  x=1:nfast;
  w1=0.5;
  w2=0.25;

  b1=bar(x, weight(2,1:nfast), w1);
  hold on;
  b2=bar(x, weight(1,1:nfast), w2);
  title(['ndtfast = ', num2str(ndtfast), blanks(5),                     ...
         'nfast = ',   num2str(nfast)]);
  grid on;
  hold off;
end

return
