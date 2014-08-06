function [ubar,vbar]=uv_barotropic(u,v,Hz,boundary);

%
% UV_BAROTROPIC:  Compute barotropic velocity components.
%
% [ubar,vbar]=uv_barotropic(u,v,Hz,boundary)
%
% This function computes vertically integrated (barotropic)
% velocity components for ROMS full grid or boundary edges.
%
% On Input:
%
%    u           U-component velocity (array or structure, U-points)
%    v           V-component velocity (array or structure, V-points)
%    Hz          Vertical level thicknesses (m) at RHO-points (3D array)
%
%    boundary    ROMS open boundary switch (1D array, OPTIONAL):
%
%                  boundary(1)        Process western  boundary
%                  boundary(2)        Process eastern  boundary
%                  boundary(3)        Process southern boundary
%                  boundary(4)        Process northern boundary
%
% On Output:
%
%    ubar        Vertically integrated U-velocity (array or structure)
%    vbar        Vertically integrated V-velocity (array or strcuture)
%
% If processing boundary data, the velocity data are structures
% with the following fields:
%
%    u.west   v.west   ubar.west   vbar.west
%    u.east   v.east   ubar.east   vbar.east
%    u.south  v.south  ubar.south  vbar.south
%    u.north  v.norht  ubar.north  vbar.north
%

% svn $Id: uv_barotropic.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Set internal parameters.

if (nargin < 4),
  DoBoundary=0;
else,
  DoBoundary=1;
end,

if (DoBoundary),
  Bnames={'west','east','south','north'};
  ufields=isfield(u,Bnames);
  vfields=isfield(v,Bnames);
  for ib=1:4,
    if (boundary(ib)),
      if (~ufields(ib)),
        error([ 'UV_BAROTROPIC - cannot find field: u.',char(Bnames(ib))]);
      end,
      if (~vfields(ib)),
        error([ 'UV_BAROTROPIC - cannot find field: v.',char(Bnames(ib))]);
      end,
    end,
  end,
end,

[Lr,Mr,N]=size(Hz);

Lu=Lr-1;  Mu=Mr;
Lv=Lr;    Mv=Mr-1;

%----------------------------------------------------------------------------
%  Vertically integrate velocity components.
%----------------------------------------------------------------------------

%  Processing grid boundary edges.

if (DoBoundary),

  for ib=1:4,
    if (boundary(ib)),
      switch ib
        case {1}
          Du=squeeze(zeros(size(u.west(:,1))));
          Dv=squeeze(zeros(size(v.west(:,1))));

          ubar.west=zeros(size(Du));
          vbar.west=zeros(size(Dv));
          for k=1:N,
            Duk=squeeze(0.5.*(Hz(1,:,k)+Hz(2,:,k)))';
            Du=Du+Duk;
            ubar.west=ubar.west+Duk.*squeeze(u.west(:,k));

            Dvk=squeeze(0.5.*(Hz(1,1:Mr-1,k)+Hz(1,2:Mr,k)))';
            Dv=Dv+Dvk;
            vbar.west=vbar.west+Dvk.*v.west(:,k);
          end,
          ubar.west=ubar.west./Du;
          vbar.west=vbar.west./Dv;
        case {2}
          Du=squeeze(zeros(size(u.east(:,1))));
          Dv=squeeze(zeros(size(v.east(:,1))));

          ubar.east=zeros(size(Du));
          vbar.east=zeros(size(Dv));
          for k=1:N,
            Duk=squeeze(0.5.*(Hz(end-1,:,k)+Hz(end,:,k)))';
            Du=Du+Duk;
            ubar.east=ubar.east+Duk.*u.east(:,k);

            Dvk=squeeze(0.5.*(Hz(end,1:Mr-1,k)+Hz(end,2:Mr,k)))';
            Dv=Dv+Dvk;
            vbar.east=vbar.east+Dvk.*v.east(:,k);
          end,
          ubar.east=ubar.east./Du;
          vbar.east=vbar.east./Dv;
        case {3}
          Du=squeeze(zeros(size(u.south(:,1))));
          Dv=squeeze(zeros(size(v.south(:,1))));

          ubar.south=zeros(size(Du));
          vbar.south=zeros(size(Dv));
          for k=1:N,
            Duk=squeeze(0.5.*(Hz(1:Lr-1,1,k)+Hz(2:Lr,1,k)));
            Du=Du+Duk;
            ubar.south=ubar.south+Duk.*u.south(:,k);

            Dvk=squeeze(0.5.*(Hz(:,1,k)+Hz(:,2,k)));
            Dv=Dv+Dvk;
            vbar.south=vbar.south+Dvk.*v.south(:,k);
          end,
          ubar.south=ubar.south./Du;
          vbar.south=vbar.south./Dv;
        case {4}
          Du=squeeze(zeros(size(u.north(:,1))));
          Dv=squeeze(zeros(size(v.north(:,1))));

          ubar.north=zeros(size(Du));
          vbar.north=zeros(size(Dv));
          for k=1:N,
            Duk=squeeze(0.5.*(Hz(1:Lr-1,end,k)+Hz(2:Lr,end,k)));
            Du=Du+Duk;
            ubar.north=ubar.north+Duk.*u.north(:,k);

            Dvk=squeeze(0.5.*(Hz(:,end-1,k)+Hz(:,end,k)));
            Dv=Dv+Dvk;
            vbar.north=vbar.north+Dvk.*v.north(:,k);
          end,
          ubar.north=ubar.north./Du;
          vbar.north=vbar.north./Dv;
      end,
    end,
  end,

%  Processing full grid.

else,

  Du=zeros(size(u(:,:,1)));
  Dv=zeros(size(v(:,:,1)));

  ubar=zeros(size(u(:,:,1)));
  vbar=zeros(size(v(:,:,1)));
  for k=1:N,
    Duk=0.5.*(Hz(1:Lr-1,1:Mr,k)+Hz(2:Lr,1:Mr,k));
    Du=Du+Duk;
    ubar=ubar+Duk.*u(:,:,k);

    Dvk=0.5.*(Hz(1:Lr,1:Mr-1,k)+Hz(1:Lr,2:Mr,k));
    Dv=Dv+Dvk;
    vbar=vbar+Dvk.*v(:,:,k);
  end,
  ubar=ubar./Du;
  vbar=vbar./Dv;

end,

return
  

