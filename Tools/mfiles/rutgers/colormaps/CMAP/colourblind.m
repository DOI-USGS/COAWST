% COLOURBLIND - Simulate colour appearance for colour blind viewers
%
% Usage: rgbblind = colourblind(rgb, pdt)
%
% Arguments:
%       rgb - RGB image or colour map.  Must be floating point in the range
%             0-1 or UINT8 in the range  0-255.
%       pdt - String 'protanopic' or 'deuteranopic' or 'tritanopic' indicating
%             type of colour blindness to simulate.
%
% Returns:
%  rgbblind - RGB colour image with floating point values in the range 0-1
%             that simulates the specified form of colour blindness.
%
% This is an implementation of:
% Computerized simulation of color appearance for dichromats
% Hans Brettel, Franc,oise Vie'not  and John D. Mollon
% Vol. 14, No. 10/October 1997/J. Opt. Soc. Am. A  pp 2647-2655
%
% See also: COLOURBLINDLABSPACE, COLOURBLINDLMSSPACE, XYZ2LMS, LMS2RGB

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

function rgbblind = colourblind(rgb, pdt)
    
    if isa(rgb, 'uint8')  % Convert uint8 to double and rescale
        rgb = double(rgb)/255;
    end
    
    if max(rgb(:)) > 1.001
        error('Image data must be double 0-1, or uint8 0-255');
    end
        
    % Define key colours that are perceived in a common way.
    % 575nm is perceived as same yellow by trichromats, protanopes and deuteranopes
    % 475nm is perceived as same blue by trichromats, protanopes and deuteranopes
    % 660nm is perceived as same red by trichomats and tritanopes    
    % 485nm is perceived as same blue-green by trichomats and tritanopes   
    % See:
    % Computerized simulation of color appearance for dichromats
    % Hans Brettel, Franc,oise Vie'not  and John D. Mollon
    % Vol. 14, No. 10/October 1997/J. Opt. Soc. Am. A  pp 2647-2655    
    
    % Values of XYZ for these wavelengths obtained from 
    % Colour and Vision Research Laboratory at UCL
    % http://cvrl.ioo.ucl.ac.uk/
    
    % CIE 1931 2-deg XYZ coordinates of key points
    XYZ475 = [0.1421    0.1126    1.0419];
    XYZ485 = [0.0580    0.1693    0.6162];    
    XYZ575 = [0.8425    0.9154    0.0018];
    XYZ660 = [0.1649    0.0610         0];
    
    % Convert to LMS cone responses
    LMS475 = xyz2lms(XYZ475);
    LMS485 = xyz2lms(XYZ485);
    LMS575 = xyz2lms(XYZ575);
    LMS660 = xyz2lms(XYZ660);
    
    % Define a point, E, in LMS space that corresponds to a maximal equal energy
    % stimulus.  I am not quite sure how to interpret what this means from the
    % paper but I am taking it as a point with equal LMS coords = [1,1,1]. This
    % also corresponds to XYZ = [1,1,1]
    E = [1,1,1];
%    E = lab2lms([100,0,0]);    % = [1.0022, 1.0241, 0.8277]

    % The origin, E and one of the key common colours are used to form a plane
    % into which an input colour is projected according to the type of colour
    % blindness.  For each type of colour blindness colour space is
    % modelled/represented/estimated by two planar regions on each side of the
    % line from the origin to E. Each half plane extends out to one of the key
    % common colours depending on the type of colour blindness.
    
    % Transfrom input colour to LMS
    Q = rgb2lms(rgb);
    
    if ndims(rgb) == 3  % Colour image
        rgbblind = zeros(size(rgb));
        [rows,cols,chan] = size(rgb);
        
        for col = 1:cols
            for row = 1:rows
                
                % Determine the point, A, in LMS space that is used to form the
                % projection plane. Note given my specification of E=[1,1,1] the
                % ratio on the right hand side is always 1
                switch pdt
                 case 'protanopic'
                  if Q(row,col,3)/Q(row,col,2) < E(3)/E(2) 
                      A = LMS575;
                  else
                      A = LMS475;
                  end
                  
                 case 'deuteranopic'
                  if Q(row,col,3)/Q(row,col,1) < E(3)/E(1) 
                      A = LMS575;
                  else
                      A = LMS475;
                  end        
                  
                 case 'tritanopic'
                  if Q(row,col,2)/Q(row,col,1) < E(2)/E(1) 
                      A = LMS660;
                  else
                      A = LMS485;
                  end        
                  
                 otherwise
                  error('Undefined colour blindness specification')
                end
                
                % Determine coefficients of the plane in LMS space
                % a*L + b*M +c*S = 0
                a = E(2)*A(3) - E(3)*A(2);
                b = E(3)*A(1) - E(1)*A(3);
                c = E(1)*A(2) - E(2)*A(1);
                
                % Perform projection of the colour in LMS space according to the
                % type of colour blindness.
                switch pdt
                 case 'protanopic'
                  Q(row,col,1) = -(b*Q(row,col,2) + c*Q(row,col,3))/a;
                 case 'deuteranopic'
                  Q(row,col,2) = -(a*Q(row,col,1) + c*Q(row,col,3))/b;
                 case 'tritanopic'
                  Q(row,col,3) = -(a*Q(row,col,1) + b*Q(row,col,2))/c;
                end
                
            end
        end
        
    elseif ndims(rgb) == 2 && size(rgb,2) == 3 % We have a colour map
        
        rgbblind = zeros(size(rgb));
        ncolours = size(rgb,1);
        
        for n = 1:ncolours
            switch pdt
             case 'protanopic'
              if Q(n,3)/Q(n,2) < E(3)/E(2) 
                  A = LMS575;
              else
                  A = LMS475;
              end
              
             case 'deuteranopic'
              if Q(n,3)/Q(n,1) < E(3)/E(1) 
                  A = LMS575;
              else
                  A = LMS475;
              end        
              
             case 'tritanopic'
              if Q(n,2)/Q(n,1) < E(2)/E(1) 
                  A = LMS660;
              else
                  A = LMS485;
              end        
              
             otherwise
              error('Undefined colour blindness specification')
            end
            
            % Determine coefficients of the plane in LMS space
            % a*L + b*M +c*S = 0
            a = E(2)*A(3) - E(3)*A(2);
            b = E(3)*A(1) - E(1)*A(3);
            c = E(1)*A(2) - E(2)*A(1);
            
            % Perform projection of the colour in LMS space according to the
            % type of colour blindness.
            switch pdt
             case 'protanopic'
              Q(n,1) = -(b*Q(n,2) + c*Q(n,3))/a;
             case 'deuteranopic'
              Q(n,2) = -(a*Q(n,1) + c*Q(n,3))/b;
             case 'tritanopic'
              Q(n,3) = -(a*Q(n,1) + b*Q(n,2))/c;
            end
            
        end
    end

    % Finally convert point Q back to RGB
    rgbblind = lms2rgb(Q);
    
    % ** Should check Q before and after projection and check for gamut
    % clipping on conversion to rgb