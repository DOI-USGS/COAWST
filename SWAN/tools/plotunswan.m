function plotunswan(matfile,basename,wavepar);
% Plots a wave parameter on unstructured grid
%
% Example:
%
% SWAN generated a binary Matlab file called 'f32har01.mat'
% TRIANGLE generated files with basename 'f32hari', e.g. 'f32hari.ele'
% To make a plot of the significant wave height, give the following
% command in Matlab:
%
%     plotunswan('f32har01','f32hari','Hsig')
%
% For other wave parameters, type the following command:
%
%    who -file f32har01
%
%
% Author  : Marcel Zijlema
% Date    : February 13, 2008
% Version : 1.0

if nargin~=3
   error('Wrong number of arguments. See "help plotunswan"')
end

eval(['load ' matfile]);                 % load binary file containing SWAN results
                                         % obtained using BLOCK command with COMPGRID-set
elefile=[basename '.ele'];
fid = fopen(elefile);                    % load TRIANGLE element based connectivity file
[nelem] = fscanf(fid,'%i',[1 3]);        % get number of triangles
ncol = 4+nelem(3);                       % specify number of columns in elefile
tri = fscanf(fid,'%i',[ncol nelem(1)])'; % get connectivity table
z=eval([wavepar]);                       % get wave parameter
trisurf(tri(:,2:4),Xp,Yp,z)              % make plot using trisurf
view(0,90);shading interp;               % make 2D view and smooth plot
colormap(jet);colorbar;axis equal        % include colorbar and equal axes
