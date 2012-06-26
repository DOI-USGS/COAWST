function [Uscope,Vscope]=uv_scope(Rscope);

%
% UV_SCOPE:  Computes adjoint sensitivity scope U- and V-mask
%
%  [Uscope,Vscope]=uv_scope(Rscope)
%
%  This function computes the adjoint sensitivity scope masks on
%  U- and V-points from the scope on RHO-points.
%
%  On Input:
%
%    Rscope        Scope mask on RHO-points (real matrix).
%
%  On Output:
%
%    Uscope        Scope mask on U-points (real matrix).
%    Vscope        Scope mask on V-points (real matrix).
%

% svn $Id: uv_scope.m 436 2010-01-02 17:07:55Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

[Lp,Mp]=size(Rscope);
L=Lp-1;
M=Mp-1;

%  Scope mask on U-points.

Uscope(1:L,1:Mp)=Rscope(2:Lp,1:Mp).*Rscope(1:L,1:Mp);

%  Scope mask on V-points.

Vscope(1:Lp,1:M)=Rscope(1:Lp,2:Mp).*Rscope(1:Lp,1:M);

return
