function [r]=rfactor(h,rmask);
 
%
% RFACTOR:  Compute bathymetry r-factor
%
% [r]=rfactor(h,rmask)
%
% This function computes the bathymetry stiffness ratio, r-factor
%
% On Input:
%
%    h           bathymetry at RHO-points (2D array)
%    rmask       Land/Sea masking at RHO-points (2D array, OPTIONAL)
%
% On Output:
%
%    r           R-factor (2D array)
%

% svn $Id: rfactor.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%


[Lp,Mp]=size(h);
L=Lp-1;
M=Mp-1;

if (nargin < 2),
  rmask=ones(size(h));
end,

%  Land/Sea mask on U-points.

for j=1:Mp,
  for i=2:Lp,
    umask(i-1,j)=rmask(i,j)*rmask(i-1,j);
  end,
end,

%  Land/Sea mask on V-points.

for j=2:Mp,
  for i=1:Lp,
    vmask(i,j-1)=rmask(i,j)*rmask(i,j-1);
  end,
end,

%----------------------------------------------------------------------------
%  Compute R-factor.
%----------------------------------------------------------------------------

hx(1:L,1:Mp)=abs(h(2:Lp,1:Mp)-h(1:L,1:Mp))./(h(2:Lp,1:Mp)+h(1:L,1:Mp));
hy(1:Lp,1:M)=abs(h(1:Lp,2:Mp)-h(1:Lp,1:M))./(h(1:Lp,2:Mp)+h(1:Lp,1:M));

hx=hx.*umask;
hy=hy.*vmask;

r(1:L,1:M)=max(max(hx(1:L,1:M),hx(1:L,2:Mp)), ...
               max(hy(1:L,1:M),hy(2:Lp,1:M)));

rmin=min(min(r));
rmax=max(max(r));
ravg=mean(mean(r));
rmed=median(median(r));

disp('  ');
disp(['Minimum r-value = ', num2str(rmin)]);
disp(['Maximum r-value = ', num2str(rmax)]);
disp(['Mean    r-value = ', num2str(ravg)]);
disp(['Median  r-value = ', num2str(rmed)]);

return
