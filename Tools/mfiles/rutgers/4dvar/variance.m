function [Fvar]=variance(Fname,Vname,Favg,Tstr,Tend);

%
% VARIANCE:  Computes the variance of requested NetCDF variable
%
% [Fvar]=variance(Fname,Vname,Favg)
%
% This function computes the variance of requested NetCDF variable from
% its specified time mean.
%
% On Input:
%
%    Fname       NetCDF file name (character string)
%    Vname       NetCDF variable name to process (character string)
%    Favg        Variable time mean (array)
%    Tstr        Starting time record to process (integer, OPTIONAL)
%    Tend        Ending   time record to process (integer, OPTIONAL)
%
% On Output:
%
%    Fvar        Requested variable variance (squared field units, array)
%

% svn $Id: variance.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Inquire number of time records.

Info = nc_vinfo(Fname,Vname);

nvdims = length(Info.Dimensions);

for n=1:nvdims,
  name = char(Info.Dimensions(n).Name);
  switch name
    case 'time',
      Nrec = dsizes(n);
  end
end

if (nargin < 4),
  Tstr = 1;
  Tend = Nrec;
end

% Read in field and compute the variance from its time mean (unbiased
% estimate since we are dividing by (ic-1).

Fvar = zeros(size(nc_read(Fname,Vname,Tstr)));

ic = 0;

for n=Tstr:Tend,
  f = nc_read(Fname,Vname,n);
  Fvar = Fvar+(f-Favg).^2;
  ic = ic+1;
end

Fvar = Fvar./max(1,ic-1);

return
