function [Fout]=shapiro1(Finp,order,scheme);

%
% SHAPIRO1:  1D Shapiro filter 
%
% [Fout]=shapiro1(Finp,order,scheme);
%
% This function smooths requested 1D array using a Shapiro filter of
% specified order.
%
% On Input:
%
%    Finp        Field be filtered (1D array)
%    order       Order of the Shapiro filter (2,4,8,16,...)
%    scheme      Switch indicating the type of boundary scheme to use:
%                  scheme = 1  =>  No change at wall, constant order
%                  scheme = 2  =>  Smoothing at wall, constant order
%                  scheme = 3  =>  No change at wall, reduced order
%                  scheme = 4  =>  Smoothing at wall, reduced order
%                  scheme = 5  =>  Periodic, constant order
%
%  On Output:
%
%     Fout       Filtered field (1D array)
%

% svn $Id: shapiro1.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

fourk=[2.500000D-1   6.250000D-2    1.562500D-2    3.906250D-3     ...
       9.765625D-4   2.44140625D-4  6.103515625D-5 1.5258789063D-5 ...
       3.814697D-6   9.536743D-7    2.384186D-7    5.960464D-8     ...
       1.490116D-8   3.725290D-9    9.313226D-10   2.328306D-10    ...
       5.820766D-11  1.455192D-11   3.637979D-12   9.094947D-13];

if (nargin < 3),
  scheme=1;
end;

Im=length(Finp);
order2=fix(order/2);

cor=zeros([1 Im]);
Fcor=zeros([1 Im]);

%----------------------------------------------------------------------------
% Compute filter correction.
%----------------------------------------------------------------------------

if (scheme == 1),

% Scheme 1:  constant order and no change at wall.

  for n=1:order2,
    if (n ~= order2),
      cor(1)=2.0*(Finp(1)-Finp(2));
      cor(Im)=2.0*(Finp(Im)-Finp(Im-1));
    else,
      cor(1)=0.0;
      cor(Im)=0.0;
    end,
    cor(2:Im-1)=2.0.*Finp(2:Im-1) - Finp(1:Im-2) - Finp(3:Im);
  end,
  Fcor=cor.*fourk(order2);

elseif (scheme == 2),

% Scheme 2:  constant order, smoothed at edges.

  for n=1:order2,
    cor(1)=2.0*(Finp(1)-Finp(2));
    cor(Im)=2.0*(Finp(Im)-Finp(Im-1));
    cor(2:Im-1)=2.0.*Finp(2:Im-1) - Finp(1:Im-2) - Finp(3:Im);
  end,
  Fcor=cor.*fourk(order2);

elseif (scheme == 3),

% Scheme 3:  reduced order and no change at wall.

  for n=1:order2,
    Istr=n;
    Iend=Im-k+1;
    if (n == 1),
      cor(2:Im-1)=2.0.*Finp(2:Im-1) - Finp(1:Im-2) - Finp(3:Im);
      cor(1)=2.0*(Finp(1)-Finp(2));
      cor(Im)=2.0*(Finp(Im)-Finp(Im-1));
    else,
      cor(Istr:Iend)=2.0.*Finp(Istr:Iend)- Finp(Istr-1:Iend-1) -  ...
                     Finp(Istr+1:Iend+1);
    end,
     Fcor(Istr)=cor(Istr)*fourk(n);
     Fcor(Iend)=cor(Iend)*fourk(n);
  end,
  Fcor(1)=0.0;
  Fcor(Istr:Iend)=cor(Istr:Iend)*fourk(order2);
  Fcor(Im)=0.0;

elseif (scheme == 4),

% Scheme 4:  reduced order, smoothed at edges.

  for n=1:order2,
    Istr=n;
    Iend=Im-k+1;
    if (n == 1),
      cor(2:Im-1)=2.0.*Finp(2:Im-1) - Finp(1:Im-2) - Finp(3:Im);
      cor(1)=2.0*(Finp(1)-Finp(2));
      cor(Im)=2.0*(Finp(Im)-Finp(Im-1));
    else,
      cor(Istr:Iend)=2.0.*Finp(Istr:Iend)- Finp(Istr-1:Iend-1) -  ...
                     Finp(Istr+1:Iend+1);
    end,
     Fcor(Istr)=cor(Istr)*fourk(n);
     Fcor(Iend)=cor(Iend)*fourk(n);
  end,
  Fcor(Istr:Iend)=cor(Istr:Iend)*fourk(order2);

elseif (scheme == 5),

% Scheme 5:  constant order, periodic.

  for n=1:order2,
    cor(2:Im-1)=2.0.*Finp(2:Im-1) - Finp(1:Im-2) - Finp(3:Im);
    cor(1)=Finp(Im-1);
    cor(Im)=Finp(2);
  end,
  Fcor=cor*fourk(order2);

end,

%----------------------------------------------------------------------------
% Apply correction.
%----------------------------------------------------------------------------

Fout=Finp-Fcor;

return
