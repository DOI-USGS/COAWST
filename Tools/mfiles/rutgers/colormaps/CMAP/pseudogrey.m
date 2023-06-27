% PSEUDOGREY Pseudogrey scale colour map with 2551 levels.
%
% Usage: map = pseudogrey
%
% Returns: map - 2551x3 array representing a pseudo-greyscale colour map with
%                2551 entries.
%
% This function implements a variation of Tyler et al's and Franzen's method for
% obtaining a pseudo-greyscale colour map with more than 256 entries.  (Their
% method produces 1786 grey levels, this one produces 2551).  This colour map
% may help if you are wanting to get the best possible rendering of a high
% dynamic range greyscale image on an 8 bit display.
%
% A typical greyscale colour map has equal RGB values increasing from 0 to 255.
%     0   0   0
%     1   1   1
%     2   2   2
%     .   .   .
%   255 255 255
%
% If we note that the conversion from RGB to grey can be approximated by
%   L = 0.3*R + 0.6*G + 0.1*B
% then by adding additional red, green or blue values to a grey colour we can
% obtain near-grey colours of intermediate greyness.  An additional red value
% changes the greyness (lightness) by approximately 0.3, green by 0.6 and blue
% by 0.1. In this implementation we obtain approximate 0.1 increments of grey as
% follows:
%
%       R G B      grey
%      _________________
%       0 0 0      0.0
%       0 0 1      0.1
%       0 0 2      0.2
%       1 0 0      0.3
%       1 0 1      0.4
%       1 0 2      0.5
%       0 1 0      0.6
%       0 1 1      0.7
%       0 1 2      0.8
%       1 1 0      0.9
%       1 1 1      1.0
%
% In practice the improvement is only just noticeable on very shallow image
% gradients.  It may possibly be useful for improving hard copy printouts on
% large posters.
%
% >> ramp = repmat(0.45:.00005:0.5, 300, 1);
% >> figure(1), imshow(ramp), colormap(gray(64))
% >> figure(2), imshow(ramp), colormap(gray(128))
% >> figure(3), imshow(ramp), colormap(gray(256))
% >> figure(4), imshow(ramp), colormap(pseudogrey)
% >> togglefigs(3,4)

% References:
%
% Christopher W. Tyler; Hoover Chan; Lei Liu; Brennan McBride; Leonid L.
% Kontsevich.  Bit stealing: how to get 1786 or more gray levels from an 8-bit
% color monitor. Proceedings Volume 1666, Human Vision, Visual Processing, and
% Digital Display III; (1992); doi: 10.1117/12.135981
%
% Also:  Rich Franzen.  PseudoGrey (1999)
% http://r0k.us/graphics/pseudoGrey.html

% Copyright (c) 2018 Peter Kovesi
% Centre for Exploration Targeting
% School of Earth Sciences
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

% PK Jan 2018

function map = pseudogrey
    
    % RGB Offsets to a grey value of equal R G and B values to obtain
    % intermediate near-grey colours at greyness increments of approximately
    % 0.1.
    offsets = [0 0 0
               0 0 1
               0 0 2
               1 0 0
               1 0 1
               1 0 2
               0 1 0
               0 1 1
               0 1 2
               1 1 0];
    
    map = zeros(2551,3);

    for n = 0:254
        map(n*10+1:(n+1)*10, 1:3) = n*ones(10,3) + offsets;
    end
               
    map(2551,:) = 255;  % Set final row
    
    % Clamp values in the last 10 rows that will have exceeded 255 due to adding
    % the offset.
    map(map>255) = 255; 
    
    map = map/255;  % Finally normalize to 0-1

