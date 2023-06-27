% XYZ2LMS - XYZ to LMS colour space conversion
%
% Usage: LMS = rgb2lms(xyz)
%
% Argument:  xyz - XYZ image or Nx3 colourmap for conversion
%
% Returns:   LMS - The converted image or colourmap.
%
% See also: RGB2LMS, LAB2RGB, RGB2NRGB, RGB2CMYK

% Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au

% PK January 2017

function LMS = xyz2lms(xyz)

    % Apply CIECAM02 transformation matrix
    % Fairchild, Mark D. (2005). Color Appearance Models (2nd Ed.). Wiley
    M = [0.7328  0.4296 -0.1624
        -0.7036  1.6975  0.0061
         0.0030  0.0136 0.9834];
    
    if ndims(xyz) == 3  % Colour image
        LMS = zeros(size(xyz));
    
        for ch = 1:3
            for n = 1:3
                LMS(:,:,ch) = LMS(:,:,ch) +  M(ch,n)*xyz(:,:,n);
            end
        end
        
    elseif size(xyz,2) == 3 % Colour map
        LMS = (M*xyz')';
    else
        error('XYZ data must be a 3 channel image or a Nx3 map');
    end
    