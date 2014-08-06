function [Fout]=shapiro2(Finp,order,scheme,napp);

%
% SHAPIRO2:  2D Shapiro filter 
%
% [Fout]=shapiro2(Finp,order,scheme);
%
% This function smooths requested 2D array using a Shapiro filter of
% specified order.
%
% On Input:
%
%    Finp        Field be filtered (2D array)
%    order       Order of the Shapiro filter (2,4,8,16,...)
%    scheme      Switch indicating the type of boundary scheme to use:
%                  scheme = 1  =>  No change at wall, constant order
%                  scheme = 2  =>  Smoothing at wall, constant order
%                  scheme = 3  =>  No change at wall, reduced order
%                  scheme = 4  =>  Smoothing at wall, reduced order
%                  scheme = 5  =>  Periodic, constant order
%    napp        Number of Shapiro filter applications (OPTIONAL)
%
% On Output:
%
%    Fout        Filtered field (2D array)
%
%  Calls:        shapiro1
%

% svn $Id: shapiro2.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

if (nargin < 4),
  napp=1;
end,

if (nargin < 3),
  scheme=1;
end;

[Im,Jm]=size(Finp);

F=Finp;

for n=1:napp,
  
%----------------------------------------------------------------------------
%  Filter all rows.
%----------------------------------------------------------------------------

  for j=1:Jm,
    Fraw=squeeze(F(:,j)); Fraw=Fraw';
    Fwrk=shapiro1(Fraw,order,scheme);
    Fout(:,j)=Fwrk';
  end,

%----------------------------------------------------------------------------
%  Filter all columns.
%----------------------------------------------------------------------------

  for i=1:Im,
    Fraw=squeeze(Fout(i,:));
    Fwrk=shapiro1(Fraw,order,scheme);
    Fout(i,:)=Fwrk;
  end,

  F=Fout;

end,
  
return
