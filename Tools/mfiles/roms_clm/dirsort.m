function D=dirsort(name)
%DIRSORT return a sorted directory listing

D=dir(name);
dn={D.name}';
[y,i]=sort(dn);
D=D(i);


