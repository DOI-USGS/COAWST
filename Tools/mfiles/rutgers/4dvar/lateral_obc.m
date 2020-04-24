function [A]=lateral_obc(A,Istr,Iend,Jstr,Jend,EW_periodic,NS_periodic);

%
% LATERAL_OBC:  Set boundary conditions for a 2D/3D field
%
%  [A]=lateral_obc(A,Istr,Iend,Jstr,Jend,EW_periodic,NS_periodic);
%
% This function applies gradient or periodic boundary conditions for
% a generic 2D/3D field.
%
% On Input:
%
%    A             Field to process (2D/3D array) 
%    Istr          Starting interior point I-index (integer)
%    Iend          Ending   interior point I-index (integer)
%    Jstr          Starting interior point J-index (integer)
%    Jend          Ending   interior point J-index (integer)
%    EW_periodic   Switch to apply E-W periodic boundary conditions
%                    (OPTIONAL)
%    NS_periodic   Switch to apply N-S periodic boundary conditions
%                    (OPTIONAL)
%
% On Output:
%
%    A             Field with updated boundary conditions (2D/3D array)
%

% svn $Id: lateral_obc.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  If not provided, turn off periodic boundary conditions.

if (nargin < 7),
  NS_periodic=0;
end,

if (nargin < 6),
  EW_periodic=0;
end,

%---------------------------------------------------------------------------
%  East-West gradient boundary conditions.
%---------------------------------------------------------------------------

if (~ EW_periodic),
  if (length(size(A)) == 2),
    A(Iend+1,Jstr:Jend)=A(Iend,Jstr:Jend);
    A(Istr-1,Jstr:Jend)=A(Istr,Jstr:Jend);
  elseif (length(size(A)) == 3),
    A(Iend+1,Jstr:Jend,:)=A(Iend,Jstr:Jend,:);
    A(Istr-1,Jstr:Jend,:)=A(Istr,Jstr:Jend,:);
  end,
end,

%---------------------------------------------------------------------------
%  North-South gradient boundary conditions.
%---------------------------------------------------------------------------

if (~ NS_periodic),
  if (length(size(A)) == 2),
    A(Istr:Iend,Jend+1)=A(Istr:Iend,Jend);
    A(Istr:Iend,Jstr-1)=A(Istr:Iend,Jstr);
  elseif (length(size(A)) == 3),
    A(Istr:Iend,Jend+1,:)=A(Istr:Iend,Jend,:);
    A(Istr:Iend,Jstr-1,:)=A(Istr:Iend,Jstr,:);
  end,
end,

%---------------------------------------------------------------------------
%  Boundary coorners.
%---------------------------------------------------------------------------

if (~ EW_periodic & ~ NS_periodic),
  if (length(size(A)) == 2),
    A(Istr-1,Jstr-1)=0.5.*(A(Istr  ,Jstr-1)+A(Istr-1,Jstr  ));
    A(Iend+1,Jstr-1)=0.5.*(A(Iend  ,Jstr-1)+A(Iend+1,Jstr  ));
    A(Istr-1,Jend+1)=0.5.*(A(Istr-1,Jend  )+A(Istr  ,Jend+1));
    A(Iend+1,Jend+1)=0.5.*(A(Iend+1,Jend  )+A(Iend  ,Jend+1));
  elseif (length(size(A)) == 3),
    A(Istr-1,Jstr-1,:)=0.5.*(A(Istr  ,Jstr-1,:)+A(Istr-1,Jstr  ,:));
    A(Iend+1,Jstr-1,:)=0.5.*(A(Iend  ,Jstr-1,:)+A(Iend+1,Jstr  ,:));
    A(Istr-1,Jend+1,:)=0.5.*(A(Istr-1,Jend  ,:)+A(Istr  ,Jend+1,:));
    A(Iend+1,Jend+1,:)=0.5.*(A(Iend+1,Jend  ,:)+A(Iend  ,Jend+1,:));
  end,
end,

%---------------------------------------------------------------------------
%  Periodic boundary conditions.
%---------------------------------------------------------------------------

if (EW_periodic | NS_periodic),
% [A]=periodic(A,Istr,Iend,Jstr,Jend,EW_periodic,NS_periodic);
end,

return
