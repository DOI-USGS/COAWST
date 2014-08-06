function [z]=depths(fname, gname, igrid, idims, tindex);

%
% DEPTHS:  Compute the ROMS depths associated with a 3D variable
%
% [z]=depths(fname, gname, igrid, idims, tindex)
%
% This function computes the depths at the requested staggered C-grid.
% If the time record (tindex) is not provided, a zero free-surface is
% assumed and the unperturbed depths are returned.
%
% On Input:
%
%    fname       NetCDF data file name (character string)
%    gname       NetCDF grid file name (character string)
%    igrid       Staggered grid C-type (integer):
%                  igrid=1  => density points
%                  igrid=2  => streamfunction points
%                  igrid=3  => u-velocity points
%                  igrid=4  => v-velocity points
%                  igrid=5  => w-velocity points
%    idims       Depths dimension order switch (integer):
%                  idims=0  => (i,j,k)  column-major order (Fortran)
%                  idims=1  => (j,i,k)  row-major order (C-language)
%    tindex      Time index (integer)
%
% On Output:
%
%    z           Depths (3D array; meters, negative)
%

% svn $Id: depths.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Check arguments and set default values.

if (nargin < 3)
  igrid=1;
end,

if (nargin < 4),
  idims=0;
end,

if (nargin < 5),
  tindex=0;
end,

% Initialize vertical tranformation and stretching function to
% original formulation values.

Vtransform=1;
Vstretching=1;
got.hc=false;

%--------------------------------------------------------------------------
% Read in S-coordinate parameters.
%--------------------------------------------------------------------------

S=nc_vnames(fname);
nvars=length(S.Variables);

for n=1:nvars
  name=char(S.Variables(n).Name);
  switch name
    case 'Vtransform'
      Vtransform=nc_read(fname,'Vtransform');
    case 'Vstretching'
      Vstretching=nc_read(fname,'Vstretching');
    case 'hc'
      got.hc=true;
    case 'sc_r'
      name_r='sc_r';
    case 's_rho'
      name_r='s_rho';
    case 'sc_w'
      name_w='sc_w';
    case 's_w'
      name_w='s_w';
  end
end

sc_r=nc_read(fname,name_r);
Cs_r=nc_read(fname,'Cs_r');

sc_w=nc_read(fname,name_w);
Cs_w=nc_read(fname,'Cs_w');

N=length(sc_r);
Np=N+1;

if (length(sc_w) == N),
  sc_w=[-1 sc_w'];
  Cs_w=[-1 Cs_w'];
end

%------------------------------------------------------------------------
% Get bottom topography.
%------------------------------------------------------------------------

h=nc_read(gname,'h');
[Lp Mp]=size(h);
L=Lp-1;
M=Mp-1;

switch ( igrid ),
  case 1
    if (idims), h=h'; end,
  case 2
    hp=0.25.*(h(1:L,1:M)+h(2:Lp,1:M)+h(1:L,2:Mp)+h(2:Lp,2:Mp));
    if (idims), hp=hp'; end,
  case 3
    hu=0.5.*(h(1:L,1:Mp)+h(2:Lp,1:Mp));
    if (idims), hu=hu'; end,
  case 4
    hv=0.5.*(h(1:Lp,1:M)+h(1:Lp,2:Mp));
    if (idims), hv=hv'; end,
  case 5
    if (idims), h=h'; end,
end

% Set critical depth parameter.

if (got.hc),
  hc=nc_read(fname,'hc');
else
  hc=min(min(h));
end

%------------------------------------------------------------------------
% Get free-surface
%------------------------------------------------------------------------

if (tindex == 0),
  zeta=zeros([Lp Mp]);
else
  zeta=nc_read(fname,'zeta',tindex);
end

switch ( igrid ),
  case 1
    if (idims), zeta=zeta'; end,
  case 2
    zetap=0.25.*(zeta(1:L,1:M )+zeta(2:Lp,1:M )+                      ...
                 zeta(1:L,2:Mp)+zeta(2:Lp,2:Mp));
    if (idims), zetap=zetap'; end,
  case 3
    zetau=0.5.*(zeta(1:L,1:Mp)+zeta(2:Lp,1:Mp));
    if (idims), zetau=zetau'; end,
  case 4
    zetav=0.5.*(zeta(1:Lp,1:M)+zeta(1:Lp,2:Mp));
    if (idims), zetav=zetav'; end,
  case 5
    if (idims), zeta=zeta'; end,
end

%------------------------------------------------------------------------
% Compute depths.
%------------------------------------------------------------------------

if (Vtransform == 1),
  switch ( igrid ),
    case 1
      for k=1:N,
	z0=(sc_r(k)-Cs_r(k))*hc + Cs_r(k).*h;
        z(:,:,k)=z0 + zeta.*(1.0 + z0./h);
      end
    case 2
      for k=1:N,
        z0=(sc_r(k)-Cs_r(k))*hc + Cs_r(k).*hp;
        z(:,:,k)=z0 + zetap.*(1.0 + z0./hp);
      end
    case 3
      for k=1:N,
        z0=(sc_r(k)-Cs_r(k))*hc + Cs_r(k).*hu;
        z(:,:,k)=z0 + zetau.*(1.0 + z0./hu);
      end
    case 4
      for k=1:N,
        z0=(sc_r(k)-Cs_r(k))*hc + Cs_r(k).*hv;
        z(:,:,k)=z0 + zetav.*(1.0 + z0./hv);
      end
    case 5
      z(:,:,1)=-h;
      for k=2:Np,
        z0=(sc_w(k)-Cs_w(k))*hc + Cs_w(k).*h;
        z(:,:,k)=z0 + zeta.*(1.0 + z0./h);
      end
  end
elseif (Vtransform == 2),
  switch ( igrid ),
    case 1
      for k=1:N,
        z0=(hc.*sc_r(k)+Cs_r(k).*h)./(hc+h);
        z(:,:,k)=zeta+(zeta+h).*z0;
      end,
    case 2
      for k=1:N,
        z0=(hc.*sc_r(k)+Cs_r(k).*hp)./(hc+hp);
        z(:,:,k)=zetap+(zetap+hp).*z0;
      end,
    case 3
      for k=1:N,
        z0=(hc.*sc_r(k)+Cs_r(k).*hu)./(hc+hu);
        z(:,:,k)=zetau+(zetau+hu).*z0;
      end,
    case 4
      for k=1:N,
        z0=(hc.*sc_r(k)+Cs_r(k).*hv)./(hc+hv);
        z(:,:,k)=zetav+(zetav+hv).*z0;
      end,
    case 5
      for k=1:Np,
        z0=(hc.*sc_w(k)+Cs_w(k).*h)./(hc+h);
        z(:,:,k)=zeta+(zeta+h).*z0;
      end
  end
end

return
