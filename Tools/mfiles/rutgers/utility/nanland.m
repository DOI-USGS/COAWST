function Fout = nanland(Finp, Ginp)

%
% NANLAND:  Mask land points to NaNs to facilitate plotting
%
% Fout = nanland(Finp, Ginp)
%
% This function mask the land point to NaNs for any ROMS array
% variable to facilitate ploting.
%
% On Input:
%
%    Finp          ROMS variable (2D or 3D array)
%
%    Ginp          ROMS Grid/History file/URL name containing all
%                    Land/Sea masking variable (string)
%
%              or, ROMS Grid structure containing all the grid
%                    variables (struct array)
%
%              or, ROMS Land/Sea masking at the appropriate C-grid
%                    location of Finp (2D array)
%
% On Output:
%
%    Fout          ROMS variable with NaN values on land points
%                    (2D or 3D array)
%
% This function can be used during plotting to mask land points.
% For example:
%
%    pcolor(rlon, rlat, nanland(salt(:,:,20),'my_grid.nc'))
% or
%    pcolor(G.lon_rho, G.lat_rho, nanland(salt(:,:,20),G))
% or
%    pcolor(G.lon_rho, G.lat_rho, nanland(salt(:,:,20),G.mask_rho))
%

% svn $Id: nanland.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Remove singleton dimension, if any.
  
Fout=squeeze(Finp);

Fsize=size(Fout);
ndims=length(Fsize);

switch ndims
  case 2
    Im=Fsize(1);
    Jm=Fsize(2);
    is3d=false;
  case 3
    Im=Fsize(1);
    Jm=Fsize(2);
    Km=Fsize(3);
    is3d=true;
end
  
%  Assign Land values to NaNs.

if (ischar(Ginp)),

  D=nc_dinfo(Ginp);
  for n=1:length(D),
    name=char(D(n).Name);
    switch name
      case 'xi_rho',
        Lr=D(n).Length;
      case 'xi_psi',
        Lp=D(n).Length;
      case 'xi_u',
        Lu=D(n).Length;
      case 'xi_v',
        Lv=D(n).Length;
      case 'eta_rho',
        Mr=D(n).Length;
      case 'eta_psi',
        Mp=D(n).Length;
      case 'eta_u',
        Mu=D(n).Length;
      case 'eta_v',
        Mv=D(n).Length;
    end
  end

  if ((Im == Lr) && (Jm == Mr)),
    mask=nc_read(Ginp,'mask_rho');
  elseif ((Im == Lp) && (Jm == Mp)),
    mask=nc_read(Ginp,'mask_psi');
  elseif ((Im == Lu) && (Jm == Mu)),
    mask=nc_read(Ginp,'mask_u');
  elseif ((Im == Lv) && (Jm == Mv)),
    mask=nc_read(Ginp,'mask_v');
  end

  if (is3d),
    ind=find(repmat(mask,[1,1,Km]) < 0.5);
  else
    ind=find(mask < 0.5);
  end

elseif (isstruct(Ginp)),

  Lr=Ginp.Lm+2; Lp=Lr-1; Lu=Lr-1; Lv=Lr; 
  Mr=Ginp.Mm+2; Mp=Mr-1; Mu=Mp;   Mv=Mp-1;

  if ((Im == Lr) && (Jm == Mr) && isfield(Ginp,'mask_rho')),
    mask=Ginp.mask_rho;
  elseif ((Im == Lp) && (Jm == Mp) && isfield(Ginp,'mask_psi')),
    mask=Ginp.mask_psi;
  elseif ((Im == Lu) && (Jm == Mu) && isfield(Ginp,'mask_u')),
    mask=Ginp.mask_u;
  elseif ((Im == Lv) && (Jm == Mv) && isfield(Ginp,'mask_v')),
    mask=Ginp.mask_v;
  else
    mask=ones(Im,Jm);
  end

  if (is3d),
    ind=find(repmat(mask,[1,1,Km]) < 0.5);
  else
    ind=find(mask < 0.5);
  end

elseif (isnumeric(Ginp)),

  Msize=size(Ginp);
  if ((Msize(1) ~= Im) || (Msize(2) ~= Jm)),
    error([' LAND - array/mask dimension mismatch: ',                 ...
           ' array = [',num2str(Im),',',num2str(Jm), '], ',           ...
           ' mask = [',num2str(Msize(1)),',',num2str(Msize(2)),']']);
  end  
  if (is3d),
    ind=find(repmat(Ginp,[1,1,Km]) < 0.5);
  else
    ind=find(Ginp < 0.5);
  end
end

%  Assign Land values to NaNs.

if (~isempty(ind)),
  Fout(ind)=NaN;
end

return