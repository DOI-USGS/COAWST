function [handle]=ptiles(NtileI, NtileJ, fname, ij_draw, varargin);

%
% PTILES:  Plot ROMS parallel horizontal tile partitions
%
% [handle]=ptiles(NtileI, NtileJ, fname, ij_draw, LineType, verbose)
%
% This function plots (overlays) parallel tile partitions in
% grid units.
%
% On Input:
%
%    NtileI      Number of parallel partitions in the I-direction
%    NtileJ      Number of parallel partitions in the J-direction
%    fname       NetCDF file name (character string)
%    ij_draw     Switch to draw in (i,j) coordinates (ij_draw=true) or
%                  in (X,Y) coordinates (ij_draw=false)
%    LineType    Line symbol and color (character string, OPTIONAL)
%    verbose     display information (default=true, OPTIONAL)
%
% On Output:
%
%    handle      plot handle
%
% calls:         tile
%

% svn $Id: ptiles.m 746 2014-12-15 23:24:27Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

LineType='k-';
verbose=true;

switch numel(varargin)
  case 1
    LineType=varargin{1};
  case 2
    LineType=varargin{1};
    verbose=varargin{2};
end

%---------------------------------------------------------------------------
% Inquire information from NetCDF file.
%---------------------------------------------------------------------------

% Get NetCDF information structure.

I=nc_inq(fname);

% Inquire about file dimensions.

for n=1:length(I.Dimensions),
  name=char(I.Dimensions(n).Name);
  switch name
    case 'xi_rho',
      Im=I.Dimensions(n).Length;
    case 'eta_rho',
      Jm=I.Dimensions(n).Length;
  end
end

% Check horizontal coordinates variables.

if (any(strcmp({I.Variables.Name}, 'spherical'))),
  spherical = nc_read(fname, 'spherical');
  if (ischar(spherical)),
    if (spherical == 'T' || spherical == 't')
      spherical = true;
    else
      spherical = false;
    end
  end
else
  spherical = true;
end
spherical=false;

for n=1:length({I.Variables.Name}),
  name=char(I.Variables(n).Name);
  switch name
    case 'lon_rho',
      X=nc_read(fname, name);
    case 'lat_rho',
      Y=nc_read(fname, name);
    case 'x_rho',
      X=nc_read(fname, name);
    case 'y_rho',
      Y=nc_read(fname, name);
  end
end

% Detemine tile partition.

Ntiles=NtileI*NtileJ-1;
Mytile=0:1:Ntiles;

[Istr,Iend,Jstr,Jend]=tile(Im-2,Jm-2,NtileI,NtileJ,Mytile,verbose);

%---------------------------------------------------------------------------
%  Draw tile boundaries.
%---------------------------------------------------------------------------

%  Draw in (I,J) coordinates.

if (ij_draw),
  hold on;

  x=1:1:Im;
  y=1:1:Jm;

  if (NtileI > 1 ),
    for i=1:NtileI-1,
      s=ones(size(y)).*Iend(i);
      handle=plot(s,y,LineType);
    end
  end

  if (NtileJ > 1 ),
    for i=1:NtileJ-1,
      j=1+(i-1)*NtileI;
      s=ones(size(x)).*Jend(j);
      handle=plot(x,s,LineType);
    end
  end

%  Draw in (X,Y) coordinates.

else

  hold on;
  
  if (NtileI > 1 ),
    for i=1:NtileI-1,
      x=squeeze(X(Iend(i),:));
      y=squeeze(Y(Iend(i),:));
      handle=plot(x,y,LineType);
    end
  end

  if (NtileJ > 1 ),
    for i=1:NtileJ-1,
      j=1+(i-1)*NtileI;
      x=squeeze(X(:,Jend(j)));
      y=squeeze(Y(:,Jend(j)));
      handle=plot(x,y,LineType);
    end
  end
  
end

return
