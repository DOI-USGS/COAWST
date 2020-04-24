function [rvor]=rvorticity(u, v, pm, pn, pmask);

% RVORTICITY:  Computes relative vorticity at PSI-points.
%
% [rvor]=rvorticity(u, v, pm, pn, pmask)
%
% This function computes 2D or 3D relative vorticity from ROMS
% fields.
%
% On Input:
%
%    u           U-momentum component (2D or 3D array, m/s)
%    v           V-momentum component (2D or 3D array, m/s)
%    pm          Reciprocal XI-grid  spacing at RHO-points (1/m)
%    pn          Reciprocal ETA-grid spacing at RHO-points (1/m)
%    pmask       Land/Sea masking at PSI-points (OPTIONAL)
%
% On Output:
%
%    rvor        Relative vorticity at PSI-points  (2D or 3D array, 1/s)
%

% svn $Id: rvorticity.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set grid paramenters.

[Lp,Mp]=size(pm);
L=Lp-1;
M=Mp-1;

if (nargin < 5),
  pmask=ones(size([L M]));
end
  
om_u=2.0./(pm(2:Lp,1:Mp)+pm(1:L,1:Mp));
on_v=2.0./(pn(1:Lp,2:Mp)+pn(1:Lp,1:M));

%  Compute relative vorticity.

if (ndims(u) == 3),
  Km=size(u,3);
  for k=1:Km
    dUde=pmask(1:L,1:M).*                                               ...
         (om_u(1:L,2:Mp).*squeeze(u(1:L,2:Mp,k))-                       ...
          om_u(1:L,1:M ).*squeeze(u(1:L,1:M ,k)));
    dVdx=pmask(1:L,1:M).*                                               ...
         (on_v(2:Lp,1:M).*squeeze(v(2:Lp,1:M,k))-                       ...
          on_v(1:L ,1:M).*squeeze(v(1:L ,1:M,k)));
    rvor(1:L,1:M,k)=0.0625.*                                            ...
                    (pm(1:L,1:M )+pm(2:Lp,1:M )+                        ...
                     pm(1:L,2:Mp)+pm(2:Lp,2:Mp)).*                      ...
                    (pn(1:L,1:M )+pn(2:Lp,1:M )+                        ...
                     pn(1:L,2:Mp)+pn(2:Lp,2:Mp)).*                      ...
                    (dVdx(1:L,1:M)-dUde(1:L,1:M));
  end
else
  dUde=pmask(1:L,1:M).*                                                 ...
       (om_u(1:L,2:Mp).*u(1:L,2:Mp)-                                    ...
        om_u(1:L,1:M ).*u(1:L,1:M ));
    dVdx=pmask(1:L,1:M).*                                               ...
         (on_v(2:Lp,1:M).*v(2:Lp,1:M)-                                  ...
          on_v(1:L ,1:M).*v(1:L ,1:M));
    rvor(1:L,1:M)=0.0625.*                                              ...
                  (pm(1:L,1:M )+pm(2:Lp,1:M )+                          ...
                   pm(1:L,2:Mp)+pm(2:Lp,2:Mp)).*                        ...
                  (pn(1:L,1:M )+pn(2:Lp,1:M )+                          ...
                   pn(1:L,2:Mp)+pn(2:Lp,2:Mp)).*                        ...
                  (dVdx(1:L,1:M)-dUde(1:L,1:M));
end

return
