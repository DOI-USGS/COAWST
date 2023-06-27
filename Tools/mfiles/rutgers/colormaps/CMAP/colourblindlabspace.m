% COLOURBLINDLABSPACE  Visualization of colour blind colour spaces in Lab
%
% This function displays a representation of the protanopic/deuteranopic or
% tritanopic colour spaces within the Lab colour space. 
%
% It is intended to support the design of colour maps that lie within the colour
% spaces of protanopes/deuteranopes or tritanopes. If a colour map lies within
% the relavant colour space then hopefully the colour map will induce the same
% visual experience for both those who are colour blind and those with normal
% vision.
%
% Usage: colourblindlabspace(pdt, Nslices, fig, colmap)
%
% Arguments:
%        pdt - String 'protanopic', 'deuteranopic' or 'tritanopic' specifying
%              the colour space to visualize. Note that the colour spaces of
%              protanopics and deuteranopics are (assumed to be) the same. 
%    Nslices - Optional number of horizontal slices of the full colour space
%              to display. Defaults to 0.
%        fig - Optional figure number to use. Defaults to 1.
%     colmap - Optional colour map or a cell array of colour maps for which
%              their colour map paths will be displayed on the plot.
%
% See also: COLOURBLINDLMSSPACE, COLOURBLIND

% Copyright (c) 2017 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% PK Auust 2017

function colourblindlabspace(pdt, Nslices, fig, colmap)
    
    if ~exist('pdt', 'var'), pdt = 'deuteranopic'; end    
    if ~exist('Nslices', 'var'), Nslices = 0; end    
    if ~exist('fig', 'var'), fig = 1; end

    fontsize = 14;
    fontweight = 'bold';
        
    %% LAB plot
    
    % Strategy.  Compute key points in colour blind space to define the two
    % planes out from the neutral axis that define the colour space.
    % Generate parametric surfaces for these two planes and their
    % corresponding colour images and map the images onto the planes.
    
    % Define key colours that are perceived in a common way.
    % 575nm is perceived as same yellow by trichromats, protanopes and deuteranopes
    % 475nm is perceived as same blue by trichromats, protanopes and deuteranopes
    % 660nm is perceived as same red by trichomats and tritanopes    
    % 485nm is perceived as same blue-green by trichomats and tritanopes          
    
    % CIE 1931 2-deg XYZ coordinates of key points
    XYZ475 = [0.1421    0.1126    1.0419];
    XYZ485 = [0.0580    0.1693    0.6162];    
    XYZ575 = [0.8425    0.9154    0.0018];
    XYZ660 = [0.1649    0.0610         0];

    % Generate key points in Lab space
    labE = [100, 0, 0];         % White
    lab475 = xyz2lab(XYZ475);   % Blue    
    lab485 = xyz2lab(XYZ485);   % Blue-green
    lab575 = xyz2lab(XYZ575);   % Yellow
    lab660 = xyz2lab(XYZ660);   % Red

    lab0 = [0, 0, 0];           % Black
     
    % Angles of planar wings towards each of the key colours
    yellowTheta  = atan2(lab575(3), lab575(2));
    blueTheta  = atan2(lab475(3), lab475(2));
    blueGreenTheta  = atan2(lab485(3), lab485(2));
    redTheta  = atan2(lab660(3), lab660(2));
 
    if strncmpi(pdt, 'protanopic', 3) || strncmpi(pdt, 'deuteranopic', 3)
        angles = [blueTheta, yellowTheta];
    elseif strncmpi(pdt, 'tritanopic', 3)
        angles = [blueGreenTheta, redTheta];
    else
        error('Colour blind type must be ''protanopic''/''deuteranopic'' or ''tritanopic');
    end
    
    % Generate meshgrid in L (Lightness) and C (chroma) in polar coords 
    [C, L] = meshgrid(0:127, 100:-1:0);
    [rows,cols] = size(L);
    wp = whitepoint('D65');
    
    figure(fig), clf
    axis([-127 127 -127 127  0 100]), hold on       
    
    % Build yellow and blue plane images
    for theta = angles
        
        lab = zeros(rows,cols,3);
        lab(:,:,1) = L;
        lab(:,:,2) = C*cos(theta);
        lab(:,:,3) = C*sin(theta);
        
        % Generate rgb values from lab
        rgb = applycform(lab, makecform('lab2srgb', 'AdaptedWhitePoint', wp));
        
        % Invert to reconstruct the lab values
        lab2 = applycform(rgb, makecform('srgb2lab', 'AdaptedWhitePoint', wp));
        
        % Where the reconstructed lab values differ from the specified values is
        % an indication that we have gone outside of the rgb gamut.  Apply a
        % mask to the rgb values accordingly
        mask = max(abs(lab-lab2),[],3);
        
        for n = 1:3
            rgb(:,:,n) = rgb(:,:,n).*(mask<2);  % tolerance of 2
        end
    
        % Generate xyz parameteric surface, this is simply the Lab coordinates
        x = lab(:,:,2);
        y = lab(:,:,3);
        z = lab(:,:,1);

        [x,y,z,rgb] = trimparametricsurface(x,y,z,rgb);        
        surface(x,y,z,rgb); shading interp; hold on
    end

    % Add horizontal colour space slices into the plot
    if Nslices > 0
        Levels = 0:100/(Nslices+1):100;
        Levels = Levels(2:end-1);  % Get rid of levels at 0 and 100
        
        for L = Levels
            rgb = generatelabslice(L);
            [x, y] = meshgrid(-127:127);
            z = ones(size(x))*L;
            
            [x,y,z,rgb] = trimparametricsurface(x,y,z,rgb);
            surface(x,y,z,rgb); shading interp; hold on
        end
    end
    
    % Generate axis tick values
    tickval = [-100 -50 0 50 100];
    tickcoords = tickval;
    ticklabels = {'-100'; '-50'; '0'; '50'; '100'};
    
    set(gca, 'xtick', tickcoords);
    set(gca, 'ytick', tickcoords);
    set(gca, 'xticklabel', ticklabels);
    set(gca, 'yticklabel', ticklabels);    
    
    ztickval = [0 20 40 60 80 100];
    zticklabels = {'0' '20' ' 40' '60' '80' '100'};
    set(gca, 'ztick', ztickval);
    set(gca, 'zticklabel', zticklabels);            
    
    set(gca,'Fontsize', fontsize);
    set(gca,'Fontweight', fontweight);
    
    % Label axes.  Note option for manual placement for a and b
    manual = 0;
    if ~manual
        xlabel('a', 'Fontsize', fontsize, 'FontWeight', fontweight);
        ylabel('b', 'Fontsize', fontsize, 'FontWeight', fontweight);
    else
        text(0, -170, 0, 'a', 'Fontsize', fontsize, 'FontWeight', ...
             fontweight);        
        text(155, 0, 0, 'b', 'Fontsize', fontsize, 'FontWeight', ...
             fontweight);
    end
    zlabel('L', 'Fontsize', fontsize, 'FontWeight', fontweight);
    
    axis vis3d
    grid on, box on, rotate3d on
    
    % Check if we have a colour map to plot
    if exist('colmap', 'var')
        if ~iscell(colmap)
            tmp = colmap;
            colmap = [];
            colmap{1} = tmp;
        end
        for n = 1:length(colmap) 
            linewidth = 2.5;
            dotcolour = [0.8, 0.8, 0.8];
            ds = 10;
            dotsize = 15;
            lab = rgb2lab(colmap{n});
            line(lab(:,2), lab(:,3), lab(:,1), 'linewidth', linewidth, 'color', ...
                 [0 0 0])
            
            plot3(lab(1:ds:end,2), lab(1:ds:end,3), lab(1:ds:end,1), '.', ...
                  'Color', dotcolour, 'MarkerSize', dotsize); 
        end
    end
    
    hold off        

%------------------------------------------------------------------------    
% Trim a parametric surface so that black areas are removed    
function [x2,y2,z2,img2] = trimparametricsurface(x,y,z,img)

    [rows,cols] = size(x);    
    x2 = x; y2 = y; z2 = z; img2 = img;

    gim = sum(img, 3);  % Sum over all colour channels
    
    for r = 1:rows
        % Find first and last non-zero elements in this row
        ind = find(gim(r,:));
        if ~isempty(ind)
            left = ind(1);
            right = ind(end);        
            
            % set all zero values in the row to be that of the first and last
            % non-zero values
            x2(r,1:left) = x(r,left);
            y2(r,1:left) = y(r,left);        
            z2(r,1:left) = z(r,left); 
            
            x2(r,right:end) = x(r,right);
            y2(r,right:end) = y(r,right);
            z2(r,right:end) = z(r,right);
            
            % Set all colours in the row to be that of the first and last
            % non-zero values.
            for ch = 1:3
                img2(r,1:left,ch) = img(r,left,ch);
                img2(r,right:end,ch) = img(r,right,ch);
            end
        end
    end
    
    % Find rows of gim that are all zero
    gimsum = sum(gim,2);
    
    ind = find(gimsum);
    if ~isempty(ind)
        top = ind(1);
        bottom = ind(end); 
    
        x2(1:top,:) = x2(top,1);
        y2(1:top,:) = y2(top,1);
        z2(1:top,:) = z2(top,1);
        
        x2(bottom:end,:) = x2(bottom,1);        
        y2(bottom:end,:) = y2(bottom,1);
        z2(bottom:end,:) = z2(bottom,1);
    end
    