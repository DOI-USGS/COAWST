function [ind]=ij2ind(a,i,j)
[m,n]=size(a);
ind=m*i-j+1;
