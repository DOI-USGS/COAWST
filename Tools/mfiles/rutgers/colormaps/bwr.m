function h = bwr(m)

% red-yellow-white-cyan-blue color map
%   BWR(M) returns an M-by-3 matrix containing the colormap.
%   BWR, by itself, is the same length as the current figure's
%   colormap. If no figure exists,  it creates one as long as MATLAB defautly does.
%-----------------------
% Syntax
%------------------------
%   %Example #1
%   ------------------------------
%   figure;
%   imagesc(peaks(150));
%   colormap(BWR(20)), colorbar
%   ------------------------------
%   Example #2
%   ------------------------------
%   figure
% bwr1(1)=subplot(1,2,1)
% imagesc(peaks(150), [0 10])
% bwr1(2)=subplot(1,2,2)
% imagesc(peaks(150), [0 10])
% colorbar(bwr1(1));colormap(bwr1(1),hot)
% colorbar(bwr1(2));colormap(bwr1(2),BWR)
%   ------------------------------
%   See also HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG,
%   COLORMAP, RGBPLOT.
%   ------------------------------
%------------------------
%Communication to:  
%------------------------
% Maosheng He
% https://orcid.org/0000-0001-6112-2499 
% Leibniz Institute of Atmospheric Physics  
% Schlossstrasse 6, 18225 Kuehlungsborn
%  hmq512@gmail.com
%  11-Feb., 2020

if nargin < 1, m = size(get(gcf,'colormap'),1); end %
x0=[0,1,3,4,5,7,8]/8;
r0=[0,0,0,1,1,1,.5];
g0=[0,0,1,1,1,0,0];
b0=r0(end:-1:1);
ColorCore=[r0; g0;b0]';
if mod(m,2)==1,
    h=interp1(x0,ColorCore,[0:1/(m-2):1]);
else
    h1=interp1(x0,ColorCore,[0:1/(m-2):0.5]);
    h2=interp1(x0,ColorCore,0.5+[0:1/(m-2):0.5]);
    h = [h1;h2];
end
