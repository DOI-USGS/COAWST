function [hout]=smooth_logbath(hinp,hmin,rmask,order,rlim,npass);

%
% SMOOTH_BATH:  Logarithmic Smoothing requested bathymetry data
%
% [hout]=smooth_bath(hinp,hmin,rmask,order,rlim)
%
% It applies a Shapiro filter to the logarithmic of the bathymetry
% data until the desired r-factor is reached.
%
% Transform input bathymetry to its logarithmic value:
%
%    hsmth = log(hinp(i,j)/hmin)
%
% Use smooth condition:
%
%    abs( hsmth(i,j) - hsmth(i+1,j) ) < log[(1+ rlim)/(1 - rlim)]
%
% Transform back smoothed bathymetry:
%
%    hout = hmin * exp(hsmth(i,j))
%
% On Input:
%
%    hinp        Input bathymetry to be smoothed (2D array; positive)
%    hmin        Minimum bathymetry value at the coast (m; positive)
%    rmask       Land/Sea masking at RHO-points (2D array)
%    order       Order of Shapiro filter (2,4,8)
%    rlim        Maximum r-factor allowed (0.35)
%    npass       Maximum number of Shapiro filter passes
%
% On Output:
%
%    hout        Smoothed bathymetry data (2D array)
%
% Calls:         rfactor, shapiro2
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                   %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.md                            Hernan G. Arango        %
%===========================================================================%

%  Set defaults.

if (nargin < 5),
  npass=50;
end,

if (nargin < 4),
  rlim=0.35;
end,

if (nargin < 3),
  order=2;
end,

%  Smooth bathymetry until desired r-factor limit is reached.

hout=log(hinp(i,j)/hmin);
r=rfactor(hout,rmask);
n=0;

while (n < npass && max(max(r)) > rlim)
  hsmth=shapiro2(hout,order,2);
  r=rfactor(hsmth,rmask);
  Im=size(r,1);
  Jm=size(r,2);
  for j=1:Jm,
    for i=1:Im
      if (r(i,j) > rlim)
        if (i < Im-2 & j < Jm-2)
          hout(i+1,j  )=hsmth(i+1,j  );
          hout(i+1,j+1)=hsmth(i+1,j+1);
          hout(i  ,j  )=hsmth(i  ,j  );
          hout(i  ,j+1)=hsmth(i  ,j+1);
        elseif (i < Im-2 & j == Jm-2)                % right edge
          hout(i  ,j  )=hsmth(i  ,j  );
          hout(i+1,j  )=hsmth(i+1,j  );
        elseif (i == Im-2 && j < Jm-2)               % top
          hout(i  ,j  )=hsmth(i  ,j  );
          hout(i  ,j+1)=hsmth(i  ,j+1);
        else                                         %  upper right corner
          hout(i  ,j  )=hsmth(i  ,j  );
        end
      end
    end
  end
  n=n+1;
end

% Transform back from log space.

hout=hmin*exp(hout);

disp(['Number of smoothing applications: ',num2str(n)]);

return


  
    
