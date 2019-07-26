function [var_rho]=v2rho_3d_mw(var_v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2000 IRD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% transfert a field at u points to a field at rho points
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,Lp,M]=size(var_v);
Mp=M+1;
Mm=M-1;
var_rho=zeros(N,Lp,M);
var_rho(:,:,2:M)=0.5*(var_v(:,:,1:Mm)+var_v(:,:,2:M));
var_rho(:,:,1)=var_rho(:,:,2);
var_rho(:,:,Mp)=var_rho(:,:,M);

return

