% COLOURBLINDLMSSPACE  Visualization of colour blind colour spaces in LMS
%
% This function displays a representation of the protanopic/deuteranopic or
% tritanopic colour spaces within the LMS colour space representing the cone
% responses of the eye.
%
% Usage: colourblindlmsspace(pdt, Slevels, fig)
%
% Arguments:
%        pdt - String 'protanopic', 'deuteranopic' or 'tritanopic' specifying
%              the colour space to visualize. Note that the colour spaces of
%              protanopics and deuteranopics are (assumed to be) the same. 
%    Slevels - Optional number of horizontal slices along the S axis of the
%              full colour space to display. Defaults to 0.
%        fig - Optional figure number to use. Defaults to 1.
%
% See also: COLOURBLINDLABSPACE, COLOURBLIND

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

function [x,y,z,rgb] = colourblindlmsspace(pdt, Slevels, fig)
    
    if ~exist('pdt', 'var'), pdt = 'deuteranopic'; end    
    if ~exist('Slevels', 'var'), Slevels = []; end    
    if ~exist('fig', 'var'), fig = 1; end

    fontsize = 14;
    fontweight = 'bold';
        
    % Strategy.  Compute key points in colour blind space to define the two
    % planes out from the neutral axis that define the colour space.
    % Generate parametric surfaces for these two planes and their
    % corresponding colour images and map the images onto the planes.
    
    % Define key colours that are perceived in a common way.
    % 575nm is perceived as same yellow by trichromats, protanopes and deuteranopes
    % 475nm is perceived as same blue by trichromats, protanopes and deuteranopes
    % 660nm is perceived as same red by trichomats and tritanopes    
    % 485nm is perceived as same blue-green by trichomats and tritanopes   
    % See:
    % Computerized simulation of color appearance for dichromats
    % Hans Brettel, Franc,oise Vie'not  and John D. Mollon
    % Vol. 14, No. 10/October 1997/J. Opt. Soc. Am. A  pp 2647-2655
    
    % CIE 1931 2-deg XYZ coordinates of key points
    XYZ475 = [0.1421    0.1126    1.0419];
    XYZ485 = [0.0580    0.1693    0.6162];    
    XYZ575 = [0.8425    0.9154    0.0018];
    XYZ660 = [0.1649    0.0610         0];

    % Generate key points in LMS space
    lmsE = xyz2lms([1, 1, 1]);  % White
    lms475 = xyz2lms(XYZ475);   % Blue    
    lms485 = xyz2lms(XYZ485);   % Blue-green
    lms575 = xyz2lms(XYZ575);   % Yellow
    lms660 = xyz2lms(XYZ660);   % Red
    lms0 = [0, 0, 0];           % Black
     
    keycol = {lmsE, lms475, lms485, lms575, lms660, lms0};

    if strncmpi(pdt, 'protanopic', 3) || strncmpi(pdt, 'deuteranopic', 3)
        colind = [2, 4];
    elseif strncmpi(pdt, 'tritanopic', 3)
        colind = [3, 5];
    elseif strncmpi(pdt, 'nothing', 2)        
        colind = [];
    else
        error('Colour blind type must be ''protanopic''/''deuteranopic'' or ''tritanopic');
    end
    
    % Generate meshgrid 0-1, 0-1 to parameterise each wing. 0-1 on the black to
    % E axis and 0-1 from the black-E axis out towards the colour point of
    % interest.
    [C, L] = meshgrid(0:.01:1);
    [rows,cols] = size(L);
    wp = whitepoint('D65');
    
    figure(fig), clf
    axis([0 1 0 1 0 1]), hold on       
    
    
    % Build the two plane images that come out from the neutral axis
    for col = colind

        % Generate unit vector perpendicular to E towards colour of interest
        N = cross(keycol{col}, lmsE);
        colV = cross(lmsE, N); colV = colV/norm(colV);
        
        % Build lms colour image over meshgrid (not right)
        lms = zeros(rows,cols,3);
        for r = 1:rows
            for c = 1:cols
                lms(r,c, :) = (rows-r)/rows*lmsE + (c-1)/cols*colV;
            end
        end

        % Generate rgb values from LMS
        rgb = lms2rgb(lms);
        
        % Invert to reconstruct the LMS values
        lms2 = rgb2lms(rgb);
        
        % Where the reconstructed lab values differ from the specified values is
        % an indication that we have gone outside of the rgb gamut.  Apply a
        % mask to the rgb values accordingly
        mask = max(abs(lms-lms2),[],3);
        for n = 1:3
            rgb(:,:,n) = rgb(:,:,n).*(mask<.05);  % tolerance of 2
        end
        
        % Generate xyz parameteric surface, this is simply the LMS coordinates
        x = lms(:,:,1);
        y = lms(:,:,2);
        z = lms(:,:,3);

        figure(fig)
        [x,y,z,rgb] = trimparametricsurface(x,y,z,rgb);        
        surface(x,y,z,rgb); shading interp; hold on
        
        line([0 1], [0 1], [0 1], 'color', [0 0 0])
    end
    
    
    % Add horizontal colour space slices into the plot
    if ~isempty(Slevels)
        
        [l, m] = meshgrid(0:.01:1);
        [rows,cols] = size(l);
        lms = zeros(rows,cols,3);  
        lms(:,:,1) = l;
        lms(:,:,2) = m;
        for S = Slevels
            lms(:,:,3) = ones(size(l))*S;
            rgb = lms2rgb(lms);
            lms2 = rgb2lms(rgb);
            
            % Where the reconstructed lab values differ from the specified values is
            % an indication that we have gone outside of the rgb gamut.  Apply a
            % mask to the rgb values accordingly
            mask = max(abs(lms-lms2),[],3);
            for n = 1:3
                rgb(:,:,n) = rgb(:,:,n).*(mask<.02);  % tolerance of 2
            end
            x = l;
            y = m;
            z = lms(:,:,3);
            [x,y,z,rgb] = trimparametricsurface(x,y,z,rgb);
            surface(x,y,z,rgb); shading interp; hold on
        end
    end


    

    % Generate axis tick values
    tickval = [0 0.2 0.4 0.6 0.8 1.0];
    tickcoords = tickval;
    ticklabels = {'0'; '0.2'; '0.4'; '0.6'; '0.8'; '1'};
    
    set(gca, 'xtick', tickcoords);
    set(gca, 'ytick', tickcoords);
    set(gca, 'ztick', tickcoords);    
    set(gca, 'xticklabel', ticklabels);
    set(gca, 'yticklabel', ticklabels);    
    set(gca, 'zticklabel', ticklabels); 
    
    set(gca,'Fontsize', fontsize);
    set(gca,'Fontweight', fontweight);
    
    xlabel('L', 'Fontsize', fontsize, 'FontWeight', fontweight);
    ylabel('M', 'Fontsize', fontsize, 'FontWeight', fontweight);    
    zlabel('S', 'Fontsize', fontsize, 'FontWeight', fontweight);

    axis vis3d
    axis equal
    grid on, box on, rotate3d on
    axis([0 1 0 1 0 1])
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
    