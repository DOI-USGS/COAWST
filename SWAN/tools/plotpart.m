function plotpart(basename);
% Plots node numbers of Metis-generated mesh partition
%
% Example:
%
% TRIANGLE generated files with basename 'f32hari', e.g. 'f32hari.node' and 'f32hari.ele'
% SWAN generated partition file 'partit.mesh' (include TEST 50,0)
% To make a plot of the colored node numbers, give the following command in Matlab:
%
%     plotpart('f32hari')
%
%
% Author  : Marcel Zijlema
% Date    : June 15, 2023
% Version : 1.0

if nargin~=1
   error('Wrong number of arguments. See "help plotpart"')
end

load partit.mesh                              % load node numbers of mesh partition
npart=max(partit);                            % get maximum node number
nodefile=[basename '.node'];
fid = fopen(nodefile);                        % load TRIANGLE vertex based connectivity file
[nnode] = fscanf(fid,'%i',[1 4]);             % get number of nodes
ncol = 3+nnode(3)+nnode(4);                   % specify number of columns in nodefile
nod = fscanf(fid,'%f',[ncol nnode(1)]);       % get mesh vertices
x=nod(2,:); y=nod(3,:);                       % get coordinates of vertices
elefile=[basename '.ele'];
fid = fopen(elefile);                         % load TRIANGLE element based connectivity file
[nelem] = fscanf(fid,'%i',[1 3]);             % get number of triangles
ncol = 4+nelem(3);                            % specify number of columns in elefile
tri = fscanf(fid,'%i',[ncol nelem(1)]);       % get connectivity table
z=partit;                                     % get node numbers
nz=size(z,1);                                 % data size for checking
if nz~=nnode(1)
   error('Partition not consistent with mesh')
end
T=tri(2:4,:);                                 % connectivity table
trisurf(T',x,y,z)                             % make plot using trisurf
view(2);                                      % make 2D view
colormap(jet(npart));axis equal               % include colorbar and equal axes
caxis([0.5 npart+0.5]);
h=colorbar;set(h,'YTick',1:npart);
