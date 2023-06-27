function varargout = cptcmap(varargin)
%CPTCMAP Apply a .cpt file as colormap to an axis
%
% cptcmap(name);
% cptcmap(name, ax);
% cptcmap(... param, val, ...);
% [cmap, lims, ticks, bfncol, ctable] = cptcmap(...)
%
% This function creates and applies a colormap defined in a color palette
% table (.cpt file).  For a full description of the cpt file format, see
% the Generic Mapping Tools documentation (http://gmt.soest.hawaii.edu/).
% Color palette files provide more flexible colormapping than Matlab's
% default schemes, including both discrete and continuous gradients, as
% well as easier direct color mapping.
%
% Limitations: X11 color names not supported, patterns not supported, CMYK
% not supported yet
%
% Input variables:
%
%   name:       .cpt file name.  You may either specify either the full
%               path, or just the file name.  In the latter case, the
%               function will look for the file in the folder specified by
%               the cptpath variable in the first line of code; by default
%               this folder is located in the same location as cptcmap.m
%               and is called cptfiles.
%
%   ax:         handle of axis or axes where colormap should be applied
%               (for pre-2014b versions, colormaps will effect the entire
%               figure(s), but axis clim adjustments for direct scaling
%               will only affect the specified axes; for 2014b+, colormap
%               will only be applied to the specified axes).  If no axis is
%               specified and no output variables are supplied, colormap
%               will be applied to the current axis.  If no axis is
%               specified and output variables are supplied, the colormap
%               will not be applied to any axes.
%
%   'showall':  When this option is used, a figure is created displaying
%               colorbars for all colormaps contained in the .cpt folder.
%               Color limits of each colormap are listed along with the
%               names of each.  A small tick mark indicates the location of
%               0, where applicable.  NOTE: the number of columns to use
%               for display is hard-coded.  As you start collecting more
%               color palettes, the figure may get too cluttered and you
%               may have to adjust this (variable ncol in the plotcmaps
%               subfunction). 
%
% Optional input variables (passed as parameter/value pairs):
%
%   'mapping':  'scaled' or 'direct'.  Scaled mapping spreads the colormap
%               to cover the color limits of the figure.  Direct mapping
%               resets the color limits of the axes so that colors are
%               mapped to the levels specified by the .cpt file. ['scaled']
%
%   'ncol':     number of colors in final colormap. If not included or NaN,
%               this function will try to choose the fewest number of
%               blocks needed to display the colormap as accurately as
%               possible. I have arbitrarily chosen that it will not try to
%               create more than 256 colors in the final colormap when
%               using this automatic scheme.  However, you can manually set
%               ncol higher if necessary to resolve all sharp breaks and
%               gradients in the colormap.
%
%   'flip':     if true, reverse the colormap order [false]
%
% Output variables:
%
%   cmap:       ncol x 3 colormap array
%
%   lims:       1 x 2 array holding minimum and maximum values for which
%               the colormap is defined.  
%
%   ticks:      vector of tick values specifying where colors were defined
%               in the original file
%
%   bfncol:     3 x 3 colormap array specifying the colors defined for
%               background (values lower than lowest color limit),
%               foreground (values higher than highest color limit), and
%               NaN values.  These do not affect the resulting colormap,
%               but can be applied by the user to replicate the behavior
%               seen in GMT.
%
%   ctable:     n x 8 color palette table, translated to Matlab color
%               space. Column 1 holds the lower limit of each color cell,
%               columns 2-4 the RGB values corresponding to the lower
%               limit, column 5 the upper limit of the color cell, and
%               columns 6-8 the RGB values of the upper limit.  When the
%               lower and upper colors are the same, this defines a
%               solid-colored cell; when they are different, colors are
%               linearly interpolated between the endpoints.
%
% Example:
%
%   [lat, lon, z] = satbath(10);
%   pcolor(lon, lat, z);
%   shading flat;
%   cptcmap('GMT_globe', 'mapping', 'direct');
%   colorbar; 

% Copyright 2011 Kelly Kearney
% File Exchange update: 4/12/2011

%------------------------------
% Parse input
%------------------------------

% The .cpt file folder.  By default, the cptfiles folder is located in
% the same place as the cptcmap.m file.  If you change this on your
% computer, change the cptpath definition to reflect the new location.
% 
% The use of two separate folders here is just to make my life easier when
% uploading this enry to github, so I can only distribute a small sample of
% colormaps, and keep my own ones separate.  It also means users can keep
% their own collection without worrying about their colormaps being
% overwritten when downloading updates of this function.  The cptfiles
% folder is included in the distribution; cptfiles_personal is not.

cptpath = {...
    fullfile(fileparts(which('cptcmap')), 'cptfiles')
    fullfile(fileparts(which('cptcmap')), 'cptfiles_personal')
    };
if ischar(cptpath)
    cptpath = {cptpath};
end
isfound = logical(cellfun(@(x) exist(x, 'dir'), cptpath));
cptpath = cptpath(isfound);
if ~any(isfound)
    error(['You have moved the cptfiles directory. '                    ...
           'Please modify the cptpath variable in this code to point '  ...
           'to the directory where your.cpt files are stored']);
end

if nargin < 1
    error('You must provide a colormap name');
end

% Check for 'showall' option first

if nargin == 1 && strcmp(varargin{1}, 'showall')
    plotcmaps(cptpath);
    return;
end

% Name of file

[blah, blah, ext] = fileparts(varargin{1});
if isempty(ext)
    varargin{1} = [varargin{1} '.cpt'];
end

if exist(varargin{1}, 'file')   % full filename and path given
    filename = varargin{1};
else                            % only file name given
    [blah,blah,ext] = fileparts(varargin{1});
    for ii = 1:length(cptpath)
        if ~isempty(ext)            % with extension
            filename = fullfile(cptpath{ii}, varargin{1});
        else                        % without extension
            filename = fullfile(cptpath{ii}, [varargin{1} '.cpt']);   
        end
        if exist(filename, 'file')
            break
        end
    end
         
    if ~exist(filename, 'file') % Not found in any folders
        error('Specified .cpt file not found');
    end
end

% Axes to which colormap will be applied

if nargin > 1 && all(ishandle(varargin{2}(:)))
    ax = varargin{2};
    pv = varargin(3:end);
    applycmap = true;
elseif nargout == 0
    ax = gca;
    pv = varargin(2:end);
    applycmap = true;
else
    pv = varargin(2:end);
    applycmap = false;
end

% Optional paramter/value pairs
    
p = inputParser;
p.addParamValue('mapping', 'scaled',                                    ...
                @(x) any(strcmpi(x, {'scaled', 'direct'})));
p.addParamValue('ncol', NaN, @(x) isscalar(x) && isnumeric(x));
p.addParamValue('flip', false, @(x) isscalar(x) && islogical(x));

p.parse(pv{:});
Opt = p.Results;
     
%------------------------------
% Calculate colormap and apply
%------------------------------

[cmap, lims,ticks,bfncol,ctable] = cpt2cmap(filename, Opt.ncol);
if Opt.flip
    if strcmp(Opt.mapping, 'direct')
      warning(['Flipping colormap with direct mapping may lead to ',    ...
               'odd color breaks']);
    end
    cmap = flipud(cmap);
end

if applycmap
    for iax = 1:numel(ax)
        if strcmp(Opt.mapping, 'direct')
            set(ax(iax), 'clim', lims);
        end
        colormap(ax(iax), cmap);
    end
end

%------------------------------
% Output
%------------------------------

allout = {cmap, lims, ticks, bfncol, ctable};
varargout(1:nargout) = allout(1:nargout);


%------------------------------
% Subfunction: Read colormap 
% from file
%------------------------------

function [cmap, lims, ticks, bfncol, ctable] = cpt2cmap(file, ncol)

% Read file

fid = fopen(file);
txt = textscan(fid, '%s', 'delimiter', '\n');
txt = txt{1};
fclose(fid);

isheader = strncmp(txt,'#',1);
isfooter = strncmp(txt,'B',1) | strncmp(txt,'F',1) | strncmp(txt,'N',1); 

% Extract color data, ignore labels (errors if other text found)

ctabletxt = txt(~isheader & ~isfooter);
ctable = str2num(strvcat(txt(~isheader & ~isfooter)));
if isempty(ctable)
    nr = size(ctabletxt,1);
    ctable = cell(nr,1);
    for ir = 1:nr
        ctable{ir} = str2num(strvcat(regexp(ctabletxt{ir},              ...
                                     '[\d\.-]*', 'match')))';
    end
    try 
        ctable = cell2mat(ctable);
    catch
        error('Cannot parse this format .cpt file yet');
    end 
end

% Determine which color model is used (RGB, HSV, CMYK, names, patterns,
% mixed)

[nr, nc] = size(ctable);
iscolmodline = cellfun(@(x) ~isempty(x), regexp(txt, 'COLOR_MODEL'));
if any(iscolmodline)
    colmodel = regexprep(txt{iscolmodline}, 'COLOR_MODEL', '');
    colmodel = strtrim(lower(regexprep(colmodel, '[#=]', '')));
else
    if nc == 8
        colmodel = 'rgb';
    elseif nc == 10
        colmodel = 'cmyk';
    else
        error('Cannot parse this format .cpt file yet');
    end
end
%     try
%         temp = str2num(strvcat(txt(~isheader & ~isfooter)));
%         if size(temp,2) == 8
%             colmodel = 'rgb';
%         elseif size(temp,2) == 10
%             colmodel = 'cmyk';
%         else % grayscale, maybe others
%             error('Cannot parse this format .cpt file yet');
%         end
%     catch % color names, mixed formats, dash placeholders
%         error('Cannot parse this format .cpt file yet');
%     end
% end
%     

% 
% iscmod = strncmp(txt, '# COLOR_MODEL', 13);
% 
% 
% if ~any(iscmod)
%     isrgb = true;
% else
%     cmodel = strtrim(regexprep(txt{iscmod}, '# COLOR_MODEL =', ''));
%     if strcmp(cmodel, 'RGB')
%         isrgb = true;
%     elseif strcmp(cmodel, 'HSV')
%         isrgb = false;
%     else
%         error('Unknown color model: %s', cmodel);
%     end
% end

% Reformat color table into one column of colors

cpt = zeros(nr*2, 4);
cpt(1:2:end,:) = ctable(:,1:4);
cpt(2:2:end,:) = ctable(:,5:8);

% Ticks

ticks = unique(cpt(:,1));

% Choose number of colors for output

if isnan(ncol)
    
    endpoints = unique(cpt(:,1));
    
    % For gradient-ed blocks, ensure at least 4 steps between endpoints
    
    issolid = all(ctable(:,2:4) == ctable(:,6:8), 2);
    
    for ie = 1:length(issolid)
        if ~issolid(ie)
            temp = linspace(endpoints(ie), endpoints(ie+1), 11)';
            endpoints = [endpoints; temp(2:end-1)];
        end
    end
    
    endpoints = sort(endpoints);
    
    % Determine largest step size that resolves all endpoints
    
    space = diff(endpoints);
    space = unique(space);
%   space = roundn(space, -3);    % To avoid floating point issues when
                                  % converting to integers
    space = round(space*1e3)/1e3;
    
    nspace = length(space);
    if ~isscalar(space)
        
        fac = 1;
        tol = .001;
        while 1
            if all(space >= 1 & (abs(space - round(space))) < tol)
                space = round(space);
                break;
            else
                space = space * 10;
                fac = fac * 10;
            end
        end
        
        pairs = nchoosek(space, 2);
        np = size(pairs,1);
        commonsp = zeros(np,1);
        for ip = 1:np
            commonsp(ip) = gcd(pairs(ip,1), pairs(ip,2));
        end
        
        space = min(commonsp);
        space = space/fac;
    end
            
    ncol = (max(endpoints) - min(endpoints))./space;
    ncol = min(ncol, 256);
    
end

% Remove replicates and mimic sharp breaks

isrep =  [false; ~any(diff(cpt),2)];
cpt = cpt(~isrep,:);

difc = diff(cpt(:,1));
minspace = min(difc(difc > 0));
isbreak = [false; difc == 0];
cpt(isbreak,1) = cpt(isbreak,1) + .01*minspace;

% Parse background, foreground, and nan colors

footer = txt(isfooter);
bfncol = nan(3,3);
for iline = 1:length(footer)
    if strcmp(footer{iline}(1), 'B')
        bfncol(1,:) = str2num(regexprep(footer{iline}, 'B', ''));
    elseif strcmp(footer{iline}(1), 'F')
        bfncol(2,:) = str2num(regexprep(footer{iline}, 'F', ''));
    elseif strcmp(footer{iline}(1), 'N')
        bfncol(3,:) = str2num(regexprep(footer{iline}, 'N', ''));
    end
end

% Convert to Matlab-format colormap and color limits

lims = [min(cpt(:,1)) max(cpt(:,1))];
endpoints = linspace(lims(1), lims(2), ncol+1);
midpoints = (endpoints(1:end-1) + endpoints(2:end))/2;

cmap = interp1(cpt(:,1), cpt(:,2:4), midpoints);

switch colmodel
    case 'rgb'
        cmap = cmap ./ 255;
        bfncol = bfncol ./ 255;
        ctable(:,[2:4 6:8]) = ctable(:,[2:4 6:8]) ./ 255;
        
    case 'hsv'
        cmap(:,1) = cmap(:,1)./360;
        cmap = hsv2rgb(cmap);
        
        bfncol(:,1) = bfncol(:,1)./360;
        bfncol = hsv2rgb(bfncol);
        
        ctable(:,2) = ctable(:,2)./360;
        ctable(:,6) = ctable(:,6)./360;
        
        ctable(:,2:4) = hsv2rgb(ctable(:,2:4));
        ctable(:,6:8) = hsv2rgb(ctable(:,6:8));
        
    case 'cmyk'
        error('CMYK color conversion not yet supported');
end

% Rouding issues: occasionally, the above calculations lead to values just
% above 1, which colormap doesn't like at all.  This is a bit kludgy, but
% should solve those issues

isnear1 = cmap > 1 & (abs(cmap-1) < 2*eps);
cmap(isnear1) = 1;

%------------------------------
% Subfunction: Display all
% colormaps
%------------------------------

function plotcmaps(folder)

for ii = 1:length(folder)
    if ii == 1
        Files = dir(fullfile(folder{ii}, '*.cpt'));
    else
        Files = [Files; dir(fullfile(folder{ii}, '*.cpt'))];
    end
end

nfile = length(Files);
ncol = 3; 
nr = ceil(nfile/ncol);
width = (1 - .05*2)/ncol;
height = (1-.05*2)/nr;
left = .05 + (0:ncol-1)*width;
bot = .05 + (0:nr-1)*height;

[l, b] = meshgrid(left, bot);
w = width * .8;
h = height * .4;

figure('color','w');
ax = axes('position', [0 0 1 1]);
hold on;

for ifile = 1:nfile
    [cmap,blah,blah,blah,ctable] = cptcmap(Files(ifile).name);
    [x,y,c] = ctable2patch(ctable);
    
    xtick = unique(x);
    dx = max(x(:)) - min(x(:));
    
    xsc = ((x-min(xtick))./dx).*w + l(ifile);
    ysc = y.*h + b(ifile);
    
    xrect = [0 1 1 0 0] .*w + l(ifile);
    yrect = [1 1 0 0 1] .*h + b(ifile);
    
    xticksc = ((xtick-min(xtick))./dx).*w + l(ifile);
    x0 = interp1(xtick, xticksc, 0);
    y0 = b(ifile) + [0 .2*h NaN .8*h h];
    x0 = ones(size(y0))*x0;

    lbl = sprintf('%s [%g, %g]',                                        ...
                  regexprep(Files(ifile).name,'\.cpt',''),              ...
                  min(x(:)), max(x(:)));
    
    patch(xsc, ysc, c, 'edgecolor', 'none');
    line(xrect, yrect, 'color', 'k');
    line(x0, y0, 'color', 'k');
    text(l(ifile), b(ifile)+h, lbl, 'interpreter', 'none',              ...
         'fontsize', 10, 'verticalalignment', 'bottom',                 ...
         'horizontalalignment', 'left');
end

set(ax, 'ylim', [0 1], 'xlim', [0 1], 'visible', 'off');

% Determine patch coordinates

function [x,y,c] = ctable2patch(ctable)

np = size(ctable,1);

x = zeros(4, np);
y = zeros(4, np);
c = zeros(4, np, 3);

y(1:2,:) = 1;

for ip = 1:np
    x(:,ip) = [ctable(ip,1) ctable(ip,5) ctable(ip,5) ctable(ip,1)];
    c(:,ip,:) = [ctable(ip,2:4);                                        ...
                 ctable(ip,6:8);                                        ...
                 ctable(ip,6:8);                                        ...
                 ctable(ip,2:4)];
end

return
