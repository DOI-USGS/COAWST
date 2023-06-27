function h = bwr2(m)

% red-white-blue color map
%   BWR2(M) returns an M-by-3 matrix containing the colormap.
%   BWR2, by itself, is the same length as the current figure's
%   colormap. If no figure exists,  it creates one as long as MATLAB defautly does.
%-----------------------
% Syntax
%------------------------
%   %Example #1
%   ------------------------------
%   figure;
%   imagesc(peaks(150));
%   colormap(BWR2(20)), colorbar
%   ------------------------------
%   Example #2
%   ------------------------------
%   figure
% bwr1(1)=subplot(1,2,1)
% imagesc(peaks(150), [0 10])
% bwr1(2)=subplot(1,2,2)
% imagesc(peaks(150), [0 10])
% try,colormap(bwr1(1),BWR), title(bwr1(1),'BWR');end;colorbar(bwr1(1));
% colormap(bwr1(2),BWR2);colorbar(bwr1(2)); title(bwr1(2),'BWR2')
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

if nargin < 1, m = size(get(gcf,'colormap'),1); end
x0=[0,1,2,3,4]/4;
ColorCore= [0 0 0.5;0 0.5 1;1 1 1;1 0 0;0.5 0 0];
if mod(m,2)==1,
    h=interp1(x0,ColorCore,[0:1/(m-2):1]);
else
    h1=interp1(x0,ColorCore,[0:1/(m-2):0.5]);
    h2=interp1(x0,ColorCore,0.5+[0:1/(m-2):0.5]);
    h = [h1;h2];
end
