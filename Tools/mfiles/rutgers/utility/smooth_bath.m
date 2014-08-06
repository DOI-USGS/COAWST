function [hout]=smooth_bath(hinp,rmask,order,rlim,npass);

%
% SMOOTH_BATH:  Smooths requested bathymetry data
%
% [hout]=smooth_bath(hinp,rmask,order,rlim)
%
% This function applies a Shapiro filter to the bathymetry data until
% the desired r-factor is reached
%
% On Input:
%
%    hinp        Input bathymetry to be smoothed (2D array)
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

% svn $Id: smooth_bath.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
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

hout=hinp;
r=rfactor(hout,rmask);
n=0;

while (n < npass && max(max(r)) > rlim),
  hsmth=shapiro2(hout,order,2);
  r=rfactor(hsmth,rmask);
  Im=size(r,1);
  Jm=size(r,2);
  for j=1:Jm,
    for i=1:Im,
      if (r(i,j) > rlim),
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
        end,
      end,
    end,
  end,
  n=n+1;
end,

disp(['Number of smoothing applications: ',num2str(n)]);

return


  
    
