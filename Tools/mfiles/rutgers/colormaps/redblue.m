function y = redblue(varargin)
% generate a RED-BLUE colormap with zero as white.
% Syntax: y = redblue(n,clim,'black')
%   Typical usage: colormap(redblue(64))
%   Positive values are displayed as blue intensities and negative
%   values displayed as red intensities. Zero is white.
%   The clim values are used to find zero to set that to white.
% Arguments:
%   All arguments are optional and can be in any order.
%   n - number of color levels 
%       (optional with default as # of colors of current colormap)
%   clim - two element vector specifying the color limits
%       (optional with default as current axis color limits)
%   black - string ('k' or 'black') specifying zero as black.
%       (optional with default as zero is white)
% See: colormap
%
% Notes:
%   This creates a custom colormap for any image so that the value 
%   zero is either white or black. The colorbar scale will be skewed
%   toward red or blue depending on the caxis values of the image.
%
%   This version flattens the edges of the colormap to improve the 
%   visualization of the gradient. This is better for larger values of n.
%
%   Keep in mind that if the scale is very skewed, there will not
%   be much of a color gradient. The gradient can always be increased
%   by using your own clim values.
%   Example:
%   y = caxis;                      % e.g. y = [-11,-5] 
%   colormap(redblue(64)            % and not much gradient
%   colormap(redblue(64,[-11,0]))   % white is -5 with larger gradient
% created   3/24/2020   Mirko Hrovat    mihrovat@outlook.com
% modified  5/04/2020   Changed algorithm for colormap calculation, eliminating
%   an error which occurs for small "clim" differences.
ntotmax = 1e4;
n = [];     clim = [];      black = false;
for k = 1:nargin
    a = varargin{k};
    switch true
        case ischar(a)
            switch true
                case strcmpi(a,'k')
                    black = true;
                case strcmpi(a,'black')
                    black = true;
            end
        case isnumeric(a) && numel(a) == 1
            n = a;
        case isnumeric(a) && numel(a) == 2
            clim = a;
    end
end
if isempty(n)
    n = size(colormap,1);
end
if isempty(clim)
    clim = caxis;
end
cmin = min(clim);
cmax = max(clim);
% np controls the visual flattening of extreme regions.
% It is the number of color levels used in these extreme regions.
% 1/8 of the range (0 to 1) uses a different gradient slope. 
% Extreme regions are at max red, min red, max blue, min blue.
% np=n/16 provides a uniform gradient but is less appealing visually.
np = max(round(sqrt(n)/4),1);   % scales slowly with n               
switch true
    case cmin >= 0      % display just blue
        ntot = round(n*cmax/(cmax-cmin));
        if ntot > ntotmax   % check if cmax is close to cmin
            ntot = ntotmax;
        end
        n2 = ntot - 2*np;
        v = [linspace(1,.875,np),linspace(.875,.125,n2),...
            linspace(.125,0,np)];
        y = repmat(v',[1,3]);
        y(:,3) = 1;
        y(1:ntot-n,:) = [];
    case cmax <= 0      % display just red
        ntot = round(n*abs(cmin/(cmax-cmin)));
        if ntot > ntotmax   % check if cmax is close to cmin
            ntot = ntotmax;
        end
        n2 = ntot - 2*np;
        v = [linspace(0,.125,np),linspace(.125,.875,n2),...
            linspace(.875,1,np)];
        y = repmat(v',[1,3]);
        y(:,1) = 1;
        y(n+1:end,:) = [];
    otherwise           % display both red and blue
        ntot = round(n*2*max(cmax,abs(cmin))/(cmax-cmin));
        if ntot > ntotmax   % check if cmax is close to cmin
            ntot = ntotmax;
        end
        m = fix(ntot/2);
        n2 = round((ntot-4*np)/2);
        v = [linspace(-1,-.875,np),linspace(-.875,-.125,n2),...
            linspace(-.125,.125,2*np),linspace(.125,.875,n2),...
            linspace(.875,1,np)];
        y = 1-abs(repmat(v',[1,3]));
        y(end-m+1:end,3) = 1;   % set blue
        y(1:m,1) = 1;           % set red
        if cmax > abs(cmin)
            y(1:ntot-n,:) = [];
        else
            y(n+1:end,:) = [];
        end   
end
if black
    y2 = y;
    y(:,1) = (1-y2(:,2)).*(y2(:,1)==1);
    y(:,3) = (1-y2(:,2)).*(y2(:,3)==1);
    y(:,2) = 0;
end
end