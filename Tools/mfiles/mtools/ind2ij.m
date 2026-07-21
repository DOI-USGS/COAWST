function [i,j]=ind2ij(a,ind)
% IND2IJ returns I,J indices of array 
%   [i,j]=ind2ij(a,ind)
%
[m,n]=size(a);
j=ceil(ind/m);
i=rem(ind,m);
%if(i==0),i=m;end;
i(i==0)=m;
