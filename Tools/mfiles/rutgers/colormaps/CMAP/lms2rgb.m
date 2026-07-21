% LMS2RGB -  LMS to RGB colour space conversion
%
% Usage: rgb = lms2rgb(LMS)
%
% Argument:  LMS - LMS image or Nx3 colourmap for conversion
%
% Returns:   rgb - The converted image or colourmap.
%
% See also: LAB2LMS, LAB2RGB, RGB2NRGB, RGB2CMYK

% Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au

% PK January 2017

function rgb = lms2rgb(LMS)

    % Apply CIECAM02 transformation matrix to convert LMS to XYZ
    % Fairchild, Mark D. (2005). Color Appearance Models (2nd Ed.). Wiley
    % LMS = M*XYZ        
        
    M = [0.7328  0.4296 -0.1624
        -0.7036  1.6975  0.0061
         0.0030  0.0136 0.9834];
    Minv = inv(M);
    
    if ndims(LMS) == 3  % Colour image
        xyz = zeros(size(LMS));
    
        for ch = 1:3
            for n = 1:3
                xyz(:,:,ch) = xyz(:,:,ch) +  Minv(ch,n)*LMS(:,:,n);
            end
        end
        
    elseif size(LMS,2) == 3 % Colour map
        xyz = (Minv*LMS')';
    else
        error('LMS data must be a 3 channel image or a Nx3 map');
    end
    
    % Finally convert from xyz to rgb
    cform = makecform('xyz2srgb');    
    rgb = applycform(xyz, cform);
    
    