function [handle]=ptiles(NtileI, NtileJ, fname, LineType);

%
% PTILES:  Plot ROMS parallel horizontal tile partitions
%
% [handle]=ptiles(NtileI, NtileJ, fname)
%
% This function plots (overlays) parallel tile partitions in
% grid units.
%
% On Input:
%
%    NtileI      Number of parallel partitions in the I-direction
%    NtileJ      Number of parallel partitions in the J-direction
%    fname       NetCDF file name (character string)
%    LineType    Line symbol and color (character string, OPTIONAL)
%
% On Output:
%
%    handle      plot handle
%
% calls:         tile
%

% svn $Id: ptiles.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

if (nargin < 4),
  LineType='k-';
end,

%---------------------------------------------------------------------------
% Inquire information from NetCDF file.
%---------------------------------------------------------------------------

% Inquire about file dimensions.

D=nc_dinfo(fname);

for n=1:length(D),
  name=char(D(n).Name);
  switch name
    case 'xi_rho',
      Im=D(n).Name;
    case 'eta_rho',
      Jm=D(n).Name;
  end,
end,

% Detemine tile partition.

Ntiles=NtileI*NtileJ-1;
Mytile=0:1:Ntiles;

[Istr,Iend,Jstr,Jend]=tile(Im-2,Jm-2,NtileI,NtileJ,Mytile);

%---------------------------------------------------------------------------
%  Draw tile boundaries.
%---------------------------------------------------------------------------

hold on;

x=1:1:Im;
y=1:1:Jm;

if (NtileI > 1 ),
  for i=1:NtileI-1,
    s=ones(size(y)).*Iend(i);
    handle=plot(s,y,LineType);
  end,
end,  

if (NtileJ > 1 ),
  for i=1:NtileJ-1,
    j=1+(i-1)*NtileI;
    s=ones(size(x)).*Jend(j);
    handle=plot(x,s,LineType);
  end,
end,

return
