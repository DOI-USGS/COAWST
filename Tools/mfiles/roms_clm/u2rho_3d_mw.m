function [var_rho]=u2rho_3d_mw(var_u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2000 IRD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% transfert a field at u points to a field at rho points
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,L,Mp]=size(var_u);
Lp=L+1;
Lm=L-1;
var_rho=zeros(N,Lp,Mp);
var_rho(:,2:L,:)=0.5*(var_u(:,1:Lm,:)+var_u(:,2:L,:));
var_rho(:,1,:)=var_rho(:,2,:);
var_rho(:,Lp,:)=var_rho(:,L,:);

return

