% RGB2LMS - RGB to LMS colour space conversion
%
% Usage: LMS = rgb2lms(rgb)
%
% Argument:  rgb - RGB image or Nx3 colourmap for conversion
%
% Returns:   LMS - The converted image or colourmap.
%
% See also: LAB2RGB, RGB2NRGB, RGB2CMYK

% Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au

% PK January 2017

function LMS = rgb2lms(rgb)

    rgb = double(rgb);
    
    % Convert to XYZ 
    cform = makecform('srgb2xyz');    
    xyz = applycform(rgb, cform);
    
    % Apply CIECAM02 transformation matrix
    % Fairchild, Mark D. (2005). Color Appearance Models (2nd Ed.). Wiley
    M = [0.7328  0.4296 -0.1624
        -0.7036  1.6975  0.0061
         0.0030  0.0136  0.9834];
    
    if ndims(rgb) == 3  % Colour image
        LMS = zeros(size(rgb));
    
        for ch = 1:3
            for n = 1:3
                LMS(:,:,ch) = LMS(:,:,ch) +  M(ch,n)*xyz(:,:,n);
            end
        end
        
    elseif size(rgb,2) == 3 % Colour map
        LMS = (M*xyz')';
    else
        error('rgb data must be a 3 channel image or a Nx3 map');
    end
    