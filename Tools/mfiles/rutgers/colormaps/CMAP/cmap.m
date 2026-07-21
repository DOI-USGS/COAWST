% CMAP Library of perceptually uniform colour maps
%
% Usage:  1:  [map, name, desc] = cmap(I, param_name, value ...)
%         2:  cmap
%         3:  cmap(str)
%  
% Arguments for Usage 1:
%
%             I - A string label indicating the colour map to be generated or a
%                 string specifying a colour map name or attribute to search
%                 for.  Type 'cmap' with no arguments to get a full list of
%                 possible colour maps and their corresponding labels.
%
%   labels:  'L1' - 'L19'  for linear maps
%            'D1' - 'D13'  for diverging maps
%            'C1' - 'C11'   for cyclic maps
%            'R1' - 'R4'   for rainbow maps
%            'I1' - 'I3'   for isoluminant maps
%
%   labels for generating maps for the colour blind:
%            'CBL1'  - 'CBL4' Linear maps for protanopic and deuteranopic viewers.
%            'CBD1'  - 'CBD2' Diverging maps for protanopic and deuteranopic viewers.
%            'CBC1'  - 'CBC2' Cyclic maps for protanopic and deuteranopic viewers.
%            'CBTL1' - 'CBTL3' Linear maps for tritanopic viewers.
%            'CBTD1'    Diverging map for tritanopic viewers.
%            'CBTC1' - 'CBTC2' Cyclic maps for tritanopic viewers.
%
%  Some colour maps have alternate labels for convenience and readability.
%  In many cases the hue string (see below) acts as an alternate label.
%  >> map = cmap('L1')  or map = cmap('grey') will produce a linear grey map.
%  >> map = cmap('D1')  or map = cmap('BWR')  blue-white-red diverging map.
%  >> cmap  % With no arguments lists all colour maps and labels.
%
%  Possible param_name - value options:
%
%     'chromaK' - The scaling to apply to the chroma values of the colour map,
%                 0 - 1.  The default is 1 giving a fully saturated colour map
%                 as designed. However, depending on your application you may
%                 want a colour map with reduced chroma/saturation values.
%                 You can use values greater than 1 however gamut clipping is
%                 likely to occur giving rise to artifacts in the colour map. 
%           'N' - Number of values in the colour map. Defaults to 256.
%       'shift' - Fraction of the colour map length N that the colour map is
%                 to be cyclically rotated, may be negative.  (Should only be
%                 applied to cyclic colour maps!). Defaults to 0.
%     'reverse' - If set to 1 reverses the colour map. Defaults to 0.
% 'diagnostics' - If set to 1 displays various diagnostic plots. Note the
%                 diagnostic plots will be for the map _before_ any cyclic
%                 shifting or reversing is applied. Defaults to 0.
%
% The parameter name strings can be abbreviated to their first letter.
%
% Returns:
%           map - Nx3 rgb colour map
%          name - A string giving a nominal name for the colour map
%          desc - A string giving a brief description of the colour map
%
%
% Usage 2 and 3:  cmap(str)
%
% Given the large number of colour maps that this function can create this usage
% option provides some help by listing the numbers of all the colour maps with
% names containing the string 'str'.  Typically this is used to search for
% colour maps having a specified attribute: 'linear', 'diverging', 'rainbow',
% 'cyclic', or 'isoluminant' etc.  If 'str' is omitted all colour maps are
% listed.  
%
%  >> cmap                % lists all colour maps
%  >> cmap('diverging')   % lists all diverging colour maps
%
% Note the listing of colour maps can be a bit slow because each colour map has to
% be created in order to determine its full name.
%
%
% Colour Map naming convention:
%
%                    linear_kryw_5-100_c67_n256
%                      /      /    |    \    \
%  Colour Map attribute(s)   /     |     \   Number of colour map entries
%                           /      |      \
%     String indicating nominal    |      Mean chroma of colour map
%     hue sequence.                |
%                              Range of lightness values
%
% In addition, the name of the colour map may have cyclic shift information
% appended to it, it may also have a flag indicating it is reversed. 
%                                              
%              cyclic_wrwbw_90-40_c42_n256_s25_r
%                                          /    \
%                                         /   Indicates that the map is reversed.
%                                        / 
%                  Percentage of colour map length
%                  that the map has been rotated by.
%
% * Attributes may be: linear, diverging, cyclic, rainbow, or isoluminant.  A
%   colour map may have more than one attribute. For example, diverging-linear or
%   cyclic-isoluminant.
%
% * Lightness values can range from 0 to 100. For linear colour maps the two
%   lightness values indicate the first and last lightness values in the
%   map. For diverging colour maps the second value indicates the lightness value
%   of the centre point of the colour map (unless it is a diverging-linear
%   colour map). For cyclic and rainbow colour maps the two values indicate the
%   minimum and maximum lightness values. Isoluminant colour maps have only
%   one lightness value. 
%
% * The string of characters indicating the nominal hue sequence uses the following code
%      r - red      g - green      b - blue
%      c - cyan     m - magenta    y - yellow
%      o - orange   v - violet 
%      k - black    w - white      j - grey
%
%   ('j' rhymes with grey). Thus a 'heat' style colour map would be indicated by
%   the string 'kryw'. If the colour map is predominantly one colour then the
%   full name of that colour may be used. Note these codes are mainly used to
%   indicate the hues of the colour map independent of the lightness/darkness and
%   saturation of the colours.
% 
% * Mean chroma/saturation is an indication of vividness of the colour map. A
%   value of 0 corresponds to a greyscale. A value of 50 or more will indicate a
%   vivid colour map.
%
% For more information about the design of these colour maps see:
% Peter Kovesi. Good Colour Maps: How to Design Them.
% https://arxiv.org/abs/1509.03700
% See also:
% https://peterkovesi.com/projects/colourmaps
%
% See also: COLORCET, EQUALISECOLOURMAP, VIEWLABSPACE, SINERAMP, CIRCLESINERAMP, COLOURMAPPATH

% Adding your own colour maps is straightforward.
%
% 1) Colour maps are almost invariably defined via a spline path through CIELAB
%    colourspace.  Use VIEWLABSPACE to work out the positions of the spline
%    control points in CIELAB space to achieve the colour map path you desire.
%    These are stored in an array 'colpts' with columns corresponding to L a and
%    b.  If desired the path can be specified in terms of RGB space by setting
%    'colourspace' to 'RGB'.  See the ternary colour maps as an example.  Note
%    the case expression for the colour map label must be upper case.
%
% 2) Set 'splineorder' to 2 for a linear spline segments. Use 3 for a quadratic
%    b-spline.
%
% 3) If the colour map path has lightness gradient reversals set 'sigma' to a
%    value of around 5 to 7 to smooth the gradient reversal.
%
% 4) If the colour map is of very low lightness contrast, or isoluminant, set
%    the lightness, a and b colour difference weight vector W to [1 1 1].
%    See EQUALISECOLOURMAP for more details
%
% 5) Set the attribute and hue sequence strings ('attributeStr' and 'hueStr')
%    appropriately so that a colour map name can be generated.  Note that if you
%    are constructing a cyclic colour map it is important that 'attributeStr'
%    contains the word 'cyclic'.  This ensures that a periodic b-spline is used
%    and also ensures that smoothing is applied in a cyclic manner.  Setting the
%    description string is optional.
%
% 6) Run CMAP specifying the number of your new colour map with the diagnostics
%    flag set to one.  Various plots are generated allowing you to check the
%    perceptual uniformity, colour map path, and any gamut clipping of your
%    colour map.
%
% Reference:
% Peter Kovesi. Good Colour Maps: How to Design Them.
% arXiv:1509.03700 [cs.GR] 2015
% https://arxiv.org/abs/1509.03700

% Copyright (c) 2013-2020 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
%
% peterkovesi.com
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% December  2013  Original version
% March     2014  Various enhancements
% August    2014  Catalogue listing
% September 2014  Provision for specifying paths in RGB, cyclic shifting and
%                 reversing 
% October   2014  Renamed to CMAP (from LABMAPLIB) and all the failed
%                 experimental maps cleaned out and retained maps renumbered
%                 in a more coherent way
% June      2016  Minor adjustments to some maps to minimise gamut clipping,
%                 this was causing some artifacts
% March     2017  Minor adjustement to L3 and L4 heatmaps. Addition of
%                 alternate names for some maps and ensured that there were
%                 description strings for all maps.
% August    2017  Addition of colour maps for colour blind
%                 viewers. Protanopic, deuteranopic and tritanopic.
% October   2017  Slight update to 'L9' to make the blue section less
%                 dominant.
% April     2018  Added L16 black-blue-green-yellow and L17, L18,
%                 L19 decreasong lightness, increasing chroma maps.
%                 Cleared out some old experimental maps.
% November  2010  Added R4, C6 and C7, reordered the naming of some cyclic maps.

function [map, name, desc] = cmap(varargin)

    [I, N, chromaK, shift, reverse, diagnostics]  = parseinputs(varargin{:});
    
    % Default parameters for colour map construction.  
    % Individual colour map specifications may override some of these.
    colourspace = 'LAB'; % Control points specified in CIELAB space
    sigma = 0;           % Default smoothing for colour map lightness equalisation.
    splineorder = 3;     % Default order of b-spline colour map curve.
    formula = 'CIE76';   % Default colour difference formula.
    W = [1 0 0];         % Colour difference weighting for Lightness,
                         % chroma and hue (default is to only use Lightness)
    desc = '';
    name = '';
    attributeStr = '';
    hueStr = '';

    % Definitions of key angles in Lab space for defining colours in colour
    % blind colour spaces. Angles are in degrees.
    yellow575 =  92.62; % For protanopic and deuteranopic colour space.
    blue475   = -79.27; 

    red660  =   32.36;  % For tritanopic colour space.
    cyan485 = -138.73; 

    
    switch I   % Note all case expressions must be upper case.
    
    %-----------------------------------------------------------------------------        
    %% Linear series
        
     case {'L1', 'L01', 'GREY', 'GRAY'}  % Grey 0 - 100 
      desc = 'Grey scale'; 
      attributeStr = 'linear';
      hueStr = 'grey';
      colpts = [  0 0 0      
                  100 0 0];
      splineorder = 2;
      
     case {'L2', 'L02', 'REDUCEDGREY'} %  Grey 10 - 95
      desc = ['Grey scale with slightly reduced contrast to '...
              'avoid display saturation problems'];
      attributeStr = 'linear';
      hueStr = 'grey';      
      colpts = [10 0 0
                95 0 0];
      splineorder = 2;
      
      % Heat map adapted from the classic definition using straight line segments from
      % black to red to yellow to white.  Here the spline order set to 3 and
      % additional control points to get a similar colour map path but with the
      % corners at red and yellow rounded off.  The best heat map out there!
     case {'L3', 'L03', 'KRYW', 'HEAT'}
      desc = 'Black-Red-Yellow-White heat colour map';
      attributeStr = 'linear';
      hueStr = 'kryw';
      colourspace = 'RGB';
      splineorder = 3;       
      colpts = [0 0 0
                .85 0 0
                1 .15 0
                1 .85 0
                1 1 .15
                1 1 1 ];      

     case {'L4', 'L04', 'KRY', 'HEATYELLOW'}
      desc = 'Black-Red-Yellow heat colour map';
      attributeStr = 'linear';
      hueStr = 'kry';
      colourspace = 'RGB';
      splineorder = 3;       
      colpts = [0 0 0
                .85 0 0
                1 .15 0
                1  1  0];

     case {'L5', 'L05', 'KGY'}
      desc = 'Colour Map along the green edge of CIELAB space';
      attributeStr = 'linear';
      hueStr = 'kgy';
      colpts = [ 5 -9  5
                 15 -23 20
                 25 -31 31
                 35 -39 39
                 45 -47 47
                 55 -55 55
                 65 -63 63
                 75 -71 71
                 85 -79 79
                 95 -38 90]; 
      
     case {'L6', 'L06', 'KBC'}
      desc = 'Blue shades running vertically up the blue edge of CIELAB space';
      attributeStr = 'linear';
      hueStr = 'kbc';
      colpts = [ 5  31 -45
                 15 50 -66
                 25 65 -90
                 35 70 -100
                 45 45 -85
                 55  20 -70
                 65  0 -53
                 75 -22 -37
                 85 -38 -20
                 95 -25 -3]; 
      
     case {'L7', 'L07', 'BMW'}
      desc = 'Blue-Pink-Light Pink colour map';
      attributeStr = 'linear';
      hueStr = 'bmw';
      colpts = [ 5 29 -43
                 15 48 -66
                 25 64 -89
                 35 73 -100
                 45 81 -86
                 55 90 -69
                 65 83 -53
                 75 56 -36
                 85 32 -22
                 95 10 -7];
      
     case {'L8', 'L08', 'BMY'}
      desc = 'Blue-Magenta-Yellow highly saturated colour map';
      attributeStr = 'linear';
      hueStr = 'bmy';
      colpts = [10 ch2ab( 55,-58)
                20 ch2ab( 75,-58)
                30 ch2ab( 75,-40)
                40 ch2ab( 73,-20)
                50 ch2ab( 75,  0)                 
                60 ch2ab( 70, 30)
                70 ch2ab( 65, 60)
                80 ch2ab( 75, 80)
                95 ch2ab( 80, 105)];
      
     case {'L9', 'L09', 'BGYW'}   % Blue - green - yellow - white
       desc = 'Blue-green-yellow colour map';
       attributeStr = 'linear';
       hueStr = 'bgyw';       
       colpts = [20  59 -80
                 35  28 -66
                 45 -14 -29
                 60 -62  60
                 85 -10  85
                 95 -15  70
                 98   0   0];
       
     case {'L10', 'GEOGRAPHIC'}
      desc = ['A "geographical" colour map.  '...
              'Best used with relief shading'];
      attributeStr = 'linear';
      hueStr = 'gow';
      colpts = [60 ch2ab(20, 180)   % pale blue green
                65 ch2ab(30, 135)
                70 ch2ab(35, 75)
                75 ch2ab(45, 85)
                80 ch2ab(22, 90)        
                85 0 0   ];
      
     case {'L11', 'GEOGRAPHIC2'}  %  Lighter version of L10 with a bit more chroma 
      desc = ['A lighter "geographical" colour map.  '...
              'Best used with relief shading'];
      attributeStr = 'linear';
      hueStr = 'gow';
      colpts = [65 ch2ab(50, 135)   % pale blue green
                75 ch2ab(45, 75)
                80 ch2ab(45, 85)
                85 ch2ab(22, 90)        
                90 0 0   ];

      
     case {'L11A', 'GEOGRAPHIC3'}  %  Lighter version of L11
      desc = ['A lighter "geographical" colour map.  '...
              'Best used with relief shading'];
      attributeStr = 'linear';
      hueStr = 'gow';
      colpts = [70 ch2ab(50, 135)   % pale blue green
                80 ch2ab(45, 75)
                85 ch2ab(45, 85)
                90 ch2ab(22, 90)        
                95 0 0   ];
      
      
     case {'L12', 'DEPTH'}
      desc =  'A "water depth" colour map';
      attributeStr = 'linear';
      hueStr = 'blue';
      colpts = [95 0 0
                80 ch2ab(20, -95)
                70 ch2ab(25, -95)
                60 ch2ab(25, -95)
                50 ch2ab(35, -95)];

      
     % The following three colour maps are for ternary images, eg Landsat images
     % and radiometric images.  These colours form the nominal red, green and
     % blue 'basis colours' that are used to form the composite image.  They are
     % designed so that they, and their secondary colours, have nearly the same
     % lightness levels and comparable chroma.  This provides consistent feature
     % salience no matter what channel-colour assignment is made.  The colour
     % maps are specified as straight lines in RGB space.  For their derivation
     % see
     % https://arxiv.org/abs/1509.03700
     
     case {'L13', 'REDTERNARY'}
      desc = 'red colour map for ternary images';
      attributeStr = 'linear';
      hueStr = 'ternary-red';
      colourspace = 'RGB';
      colpts = [0.00 0.00 0.00
                0.90 0.17 0.00];

      splineorder = 2; 
      
     case {'L14', 'GREENTERNARY'}
      desc = 'green colour map for ternary images';
      attributeStr = 'linear';
      hueStr = 'ternary-green';
      colourspace = 'RGB';
      colpts = [0.00 0.00 0.00
                0.00 0.50 0.00];
      
      splineorder = 2; 
      
     case {'L15', 'BLUETERNARY'}
      desc = 'blue colour map for ternary images';
      attributeStr = 'linear';
      hueStr = 'ternary-blue';
      colourspace = 'RGB';
      colpts = [0.00 0.00 0.00
                0.10 0.33 1.00];
      
      splineorder = 2; 
      
     % Variation of L9 that starts at black.  Works well.
     case {'L16', 'KBGYW'}   % Black - Blue - green - yellow - white
       desc = 'Black-Blue-Green-Yellow-White colour map';
       attributeStr = 'linear';
       hueStr = 'kbgyw';       
       colpts = [10   0  0
                 20  59 -80
                 35  28 -66
                 45 -14 -29
                 60 -62 60
                 85  -10  85
                 95 -15 70
                 98  0  0];
       
       W = [1,0,0];
        
      % Linear decreasing lightness, increasing chroma to blue. 
      case {'L17', 'WORB'}
        desc = 'White-Orange-Red-Blue, decreasing lightness with increasing saturation';
        attributeStr = 'linear';
        hueStr = 'worb';
        
        % Construct a spiral of increasing chroma down through Lab space
        nsteps = 12;
        ang1 = 120; ang2 = -60;   % Linearly interpolate hue angle
        ang = ang1:(ang2 - ang1)/(nsteps-1):ang2;

        % Interpolate chroma but use a 'gamma' of 0.5 to keep the colours more saturated.
        sat1 = 0; sat2 = 80;
        sat = ([sat1:(sat2-sat1)/(nsteps-1):sat2]/sat2).^0.5 * sat2;

        l1 = 100; l2 = 25;        % Linearly interpolate lightness
        l = l1:(l2-l1)/(nsteps-1):l2;

        colpts = zeros(nsteps,3);
        for n=1:nsteps
           colpts(n,:) = [l(n), ch2ab(sat(n), ang(n))];
        end
        
        splineorder = 2;        

      case {'L18', 'WYOR'}
        desc = 'White-Yellow-Orange-Red, decreasing lightness with increasing saturation';
        attributeStr = 'linear';
        hueStr = 'wyor';
        
        % Construct a spiral of increasing chroma down through Lab space
        nsteps = 12;
        ang1 = 120; ang2 = 35; % Linearly interpolate hue angle
        ang = ang1:(ang2 - ang1)/(nsteps-1):ang2;

        % Interpolate chroma but use a 'gamma' of 0.5 to keep the colours more saturated.
        sat1 = 0; sat2 = 84;
        sat = ([sat1:(sat2-sat1)/(nsteps-1):sat2]/sat2).^0.5 * sat2;

        l1 = 100; l2 = 45;        % Linearly interpolate lightness
        l = l1:(l2-l1)/(nsteps-1):l2;

        colpts = zeros(nsteps,3);
        for n=1:nsteps
           colpts(n,:) = [l(n), ch2ab(sat(n), ang(n))];
        end
        
        splineorder = 2;        
        
      case {'L19', 'WCMR'}
        desc = 'White-Cyan-Magenta-Red, decreasing lightness with increasing saturation';
        attributeStr = 'linear';
        hueStr = 'wcmr';
        
        % Construct a spiral of increasing chroma down through Lab
        % space.  Use linear interolation for hue angle, chroma and lightness
        nsteps = 12;
        ang1 = -140; ang2 = 40; 
        ang = ang1:(ang2 - ang1)/(nsteps-1):ang2;

        sat1 = 0; sat2 = 84;
        sat = sat1:(sat2-sat1)/(nsteps-1):sat2;

        l1 = 100; l2 = 45;
        l = l1:(l2-l1)/(nsteps-1):l2;

        colpts = zeros(nsteps,3);
        for n=1:nsteps
           colpts(n,:) = [l(n), ch2ab(sat(n), ang(n))];
        end
        
        splineorder = 2;        

       

        % A map inspired by Parula but with a smoother path that does not attempt to
        % include cyan and maintains a more uniform slope upwards in LAB
        % space. It also starts at a lightness of 20 with zero chroma. It works
        % very much better than Parula!
      case {'L20', 'GOULDIAN'}
        desc = 'Black-Blue-Green-Orange-Yellow map';
        attributeStr = 'linear';
        hueStr = 'kbgoy';

        colpts = [20 0  0
                  40 55 -90 
                  55 -47 0 
                  70 -20 70
                  80  20 80
                  95 -21 92];
        sigma = 0;
        splineorder = 3;
        
        

     %-------------------------------------------------------------------------
     %% Diverging colour maps
     
     % Note that on these colour maps often we do not go to full white but use a
     % lightness value of 95. This helps avoid saturation problems on monitors.
     % A lightness smoothing sigma of 5 to 7 is used to avoid generating a false
     % feature at the white point in the middle.  Note however, this does create
     % a small perceptual contrast blind spot at the middle.
        
     case {'D1', 'D01', 'BWR', 'COOLWARM'}
      desc = ['Classic diverging blue - white - red colour map. End colours are' ...
              ' matched in lightness and chroma'];
      attributeStr = 'diverging';
      hueStr = 'bwr';
      colpts = [40  ch2ab(83,-64)
                95  0   0
                40  ch2ab(83, 39)];   
      sigma = 7;
      splineorder = 2; 

      
     case {'D1A', 'D01A', 'BWRA'}
      desc = ['Diverging blue - white - red colour map. Similar to D1 but with darker'...
              ' end point colours'];
      attributeStr = 'diverging';
      hueStr = 'bwr';
      colpts = [20  ch2ab(49,-64)
                40  ch2ab(83,-64)
                65  ch2ab(60,-64)
                95  0   0
                65  ch2ab(60, 37)
                40  ch2ab(83, 37)
                20  ch2ab(49, 37)];
      sigma = 7;
      splineorder = 2; 
 
      
    case {'D1B', 'D01B', 'BWO'}
      desc = 'blue - white - dark orange colour map';
      attributeStr = 'diverging';
      hueStr = 'bwr';
      colpts = [40  ch2ab(50,-64)
                95  0   0
                40  ch2ab(50, 77)];   
      sigma = 7;
      splineorder = 2; 
      
     case {'D2', 'D02', 'GWV'}
      desc = 'Diverging green - white - violet colour map';
      attributeStr = 'diverging';
      hueStr = 'gwv';
      colpts = [55 -50  55
                95   0   0
                55  60 -55];  
      sigma = 7;
      splineorder = 2; 
     
     case {'D3', 'D03', 'GWR'}
      desc = 'Diverging green - white - red colour map';
      attributeStr = 'diverging';
      hueStr = 'gwr';
      colpts = [55 -50 55
                95   0  0
                55  63 39];  
      sigma = 7;
      splineorder = 2;  
      
     case {'D4', 'D04', 'BKR'}
      desc = 'Diverging blue - black - red colour map';
      attributeStr = 'diverging';
      hueStr = 'bkr';
      colpts = [55 ch2ab(70, -76)
                10  0   0
                55 ch2ab(70, 35)];   
      sigma = 7;
      splineorder = 2; 
      
     case {'D5', 'D05', 'GKR'}
      desc = 'Diverging green - black - red colour map';
      attributeStr = 'diverging';
      hueStr = 'gkr';
      colpts = [60 ch2ab(80, 134)
                10  0   0
                60 ch2ab(80, 40)];   
      sigma = 7;
      splineorder = 2; 
      
     case {'D6', 'D06', 'BKY'}
      desc = 'Diverging blue - black - yellow colour map';
      attributeStr = 'diverging';
      hueStr = 'bky';
      colpts = [60 ch2ab(60, -85)
                10  0   0
                60 ch2ab(60, 85)];   
      sigma = 7;
      splineorder = 2;        

     case {'D7', 'D07', 'BJY', 'DIVBJY'}  % Linear diverging  blue - grey - yellow.  Works well
      desc = ['Linear-diverging blue - grey - yellow colour map. This kind' ...
              ' of diverging map has no perceptual dead spot at the centre'];
      attributeStr = 'diverging-linear';
      hueStr = 'bjy';
      colpts = [30 ch2ab(89, -59)
                60 0 0                 
                90 ch2ab(89, 96)];
      splineorder = 2;
      
      
     case {'D8', 'D08', 'BJR'} % Linear diverging  blue - grey - red
      desc = 'Linear-diverging blue - grey - red';
      attributeStr = 'diverging-linear';
      hueStr = 'bjr';      
      colpts = [30 ch2ab(105, -58)
                42.5 0 0                 
                55 ch2ab(105,41)];
      splineorder = 2;
      
     case {'D9', 'D09'}  % Lightened version of D1 for relief shading - Good.
      desc = ['Diverging low contrast blue - white - red colour map.  Good in' ...
              ' conjunction with relief shading'];
      attributeStr = 'diverging';
      hueStr = 'bwr';
      colpts = [55  ch2ab(73,-74)
                98  0   0
                55  ch2ab(73, 39)];   
      sigma = 2;   % Less smoothing needed for low contrast
      splineorder = 2; 
      
     case 'D10' % low contrast diverging map for when you want to use
                % relief shading 
      desc = ['Diverging low contrast cyan - white - magenta colour map.  Good in' ...
              ' conjunction with relief shading'];
      attributeStr = 'diverging';
      hueStr = 'cwm';
      colpts = [80 ch2ab(44, -135)
                100  0   0
                80 ch2ab(44, -30)];   
      sigma = 0;   % No smoothing needed for lightness range of 20
      splineorder = 2; 
      W = [1 1 1];
      
     case 'D11' % Constant lightness diverging map for when you want to use
                % relief shading  ? Perhaps lighten the grey so it is not quite
                % isoluminant ?
      desc = 'Diverging isoluminant lightblue - lightgrey - orange colour map';
      attributeStr = 'diverging-isoluminant';
      hueStr = 'cjo';
      colpts = [70 ch2ab(50, -105)
                70  0   0
                70 ch2ab(50, 45)];   
      sigma = 7;
      splineorder = 2; 
      W = [1 1 1];
      
     case 'D12' % Constant lightness diverging map for when you want to use
                % relief shading  ? Perhaps lighten the grey so it is not quite
                % isoluminant ?
      desc = 'Diverging isoluminant lightblue - lightgrey - pink colour map';
      attributeStr = 'diverging-isoluminant';
      hueStr = 'cjm';
      colpts = [75 ch2ab(46, -122)
                75  0   0
                75 ch2ab(46, -30)];   
      sigma = 7;
      splineorder = 2; 
      W = [1 1 1];
      
      case {'D13', 'BWG'}  % Pleasing Blue-White-Green diverging map
        desc = 'Diverging blue - white - green colour map';
        attributeStr = 'diverging';
        hueStr = 'bwg';
        colpts = [20  ch2ab(40, -70)
                  45  ch2ab(60, -70)
                  70  ch2ab(50, -70-35)
                  95  0 0
                  95  0 0
                  70  ch2ab(50, 138+35)
                  45  ch2ab(60, 138)
                  20  ch2ab(40, 138)];
        sigma = 5;
        splineorder = 3;
 
     %-------------------------------------------------------------------------      
     %% Cyclic colour maps
    
     case {'C1', 'C01'} % I think this is my best zigzag style cyclic map - Good!
               % Control points are placed so that lightness steps up and
               % down are equalised.  Additional intermediate points are
               % placed to try to even up the 'spread' of the key colours.
      desc = ['Cyclic: magenta - red - yellow - blue. Alows four' ...
              ' orientations/phase angles to be visulaised.'];
       attributeStr = 'cyclic';
       hueStr = 'mrybm'; 
       mag = [75 60 -40];
       yel = [75 0 77];
       blu = [35  70 -100];
       red = [35 60 48];
       colpts = [mag
                 55 70 0 
                 red
                 55 35 62
                 yel
                 50 -20 -30
                 blu
                 55 45 -67
                 mag];
       
       sigma = 7;
       splineorder = 2;  % linear path 
       

     case {'C2', 'C02', 'PHASE4'}  % A big diamond across the gamut.  Really good!  Incorporates two
                            % extra cotnrol points around blue to extend the width of that
                            % segment slightly.
      desc = ['Cyclic: magenta - yellow - green - blue. Alows four' ...
              ' orientations/phase angles to be visualised.'];
      attributeStr = 'cyclic';
      hueStr = 'mygbm';
      colpts = [62.5  83 -54
                80 20 25
                95 -20 90
                62.5 -65 62   
                42 10 -50
                30 75 -103                
                48 70 -80  
                62.5  83 -54];
      sigma = 7;
      splineorder = 2;
      
     case {'C2A', 'C02A'}  % Attempt at a brighter version of C2. Not entirely successful
      desc = '';
      attributeStr = 'cyclic';
      hueStr = 'mygbm';
      colpts = [71  75 -47   % magenta
                97 -22 94    % yellow
                71 -71 68    % green
                45 41 -88    % blue  45 41 -88
                71 75 -47]; 

      sigma = 7;
      splineorder = 2;
      
      case {'C3', 'C03'}   % white-red-black-blue-white - variation of C4
      desc = ['Cyclic: white - red - black - blue'];  
      attributeStr = 'cyclic';
      hueStr = 'wrkbw';
      colpts = [90 0 0
                50 65 56
                10 0 0
                50 31 -80
                90 0 0];
      sigma = 7;
      splineorder = 2;      
      
     case {'C4', 'C04', 'PHASE2'}   % white-red-white-blue-white Works nicely
      desc = ['Cyclic: white - red - white - blue. Good if you just want to' ...
              ' visualise +ve and -ve phase'];  
      attributeStr = 'cyclic';
      hueStr = 'wrwbw';
      colpts = [90 0 0
                40 65 56
                90 0 0
                40 31 -80
                90 0 0];
      sigma = 7;
      splineorder = 2;      

     case {'C5', 'C05', 'CYCLICGREY'}   % Cyclic greyscale  Works well
      desc = 'Cyclic: greyscale'; 
      attributeStr = 'cyclic';
      hueStr = 'grey';
      colpts = [50 0 0
                85 0 0
                15 0 0
                50 0 0];
      sigma = 7;
      splineorder = 2;
      
      % A six colour map. Yellow, Cyan and Magenta at a lightness of 90. (Magenta gets
      % washed out with little chroma at 90) Red, Green and Blue at a lightness
      % of 50. The Green is a bit dark and the blue a bit light at 50. Overall a
      % slightly strange colour circle but I think it works quite well.      
      
      case {'C6', 'C06'}  
        desc = 'Six colour cyclic with primaries and secondaries matched in lightness.'; 
        attributeStr = 'cyclic';
        hueStr = 'rygcbmr';
         
         y90 = [90 -7 90];
         m90 = [90 24 -17];
         c90 = [90 -48 -14];
         
         y85 = [85 -22 84];
         m85 = [85 35 -23];
         c85 = [85 -46 -14];
         
         y75 = [75 -16 74];
         m75 = [75 60 -40];
         c75 = [75 -39 -16];         
         
         r55 = [55 80 66];
         b55 = [55 19 -71];
         g55 = [55 -60 56];
         
         r50 = [50 78 62];
         b50 = [50 30 -81];
         g50 = [50 -54 52];
         
         colpts = [r50
                   y90
                   g50
                   c90
                   b50
                   m90
                   r50];
         
         sigma = 7;
         splineorder = 2;  
         
      case {'C7', 'C07'}  % Zig-Zag Yellow - Magenta - Cyan - Green - Yellow
                          % Colours are well balanced over the quadrants 
         desc = 'Cyclic Yellow - Magenta - Cyan - Green - Yellow ';
         attributeStr = 'cyclic';
         hueStr = 'ymcgy';
         
         yel = [85 3 86];
         y90 = [90 -20 90];  % a was -6
         c90 = [90 -48 -14];
         cyn = [85 -40 -25];
         blu = [50  30 -79];
         red = [50 79 64];   
         
         g70 = [70 -71 70];
         g65 = [65 -67 66];
         g60 = [60 -64 61];
         g55 = [55 -60 57];
         
         m70 = [70 75 -47];
         m65 = [65  90 -55];         
         m60 = [60 98 -61];
         m55 = [55 94 -71];
         
         colpts = [y90
                   m60
                   c90
                   g60
                   y90];
         
       sigma = 5;
       splineorder = 2;           
       
      
     case {'C8', 'C08'}  % Elliptical path - ok
      desc = 'Cyclic: low contrast colour ellipse';
      attributeStr = 'cyclic';
      hueStr = 'mygbm';
      ang = 112;
      colpts = [70    ch2ab(46,  ang-90)
                90    ch2ab(82,  ang)
                70    ch2ab(46,  ang+90)
                50    ch2ab(82,  ang+180)
                70    ch2ab(46,  ang-90)];
      W = [1 1 1];
      
     case {'C9', 'C09'} % Elliptical path.  Greater range of lightness values and
               % slightly more saturated colours.  Seems to work however I do
               % not find the colour sequence that attractive. This is a
               % constraint of the gamut.
      desc = 'Cyclic: colour ellipse';
      attributeStr = 'cyclic';
      hueStr = 'mybm';      
      ang = 124;
      colpts = [60    ch2ab(40,  ang-90)
                100   ch2ab(98,  ang)
                60    ch2ab(40,  ang+90)
                20    ch2ab(98,  ang+180)
                60    ch2ab(40,  ang-90)];
      W = [1 1 1];
      sigma = 7;                        

       
     case 'C10'  % Circle at 67  - sort of ok but a bit flouro
      desc = 'Cyclic: isoluminant';
      attributeStr = 'cyclic-isoluminant colour circle';
      hueStr = 'mgbm';
      chr = 42;
      ang = 124;
      colpts = [67  ch2ab(chr,  ang-90)
                67  ch2ab(chr,  ang)
                67  ch2ab(chr,  ang+90)
                67  ch2ab(chr,  ang+180)
                67  ch2ab(chr,  ang-90)];
      W = [1 1 1];
      

     case 'C11' % Variation of C1. Perceptually this is good. Excellent balance
                % of colours in the quadrants but the colour mix is not to my
                % taste.  Don't like the green.  The red-green transition
                % clashes
      desc = ['Cyclic: blue - green - red - magenta. Allows four' ...
              ' orientations/phase angles to be visulaised.']; 
       attributeStr = 'cyclic';
       hueStr = 'bgrmb';       
       blu = [35  70 -100];
       colpts = [blu
                 70 -70 64
                 35 65 50
                 70 75 -46
                 blu        ];      
       sigma = 7;
       splineorder = 2;  % linear path 
         
         
     %-----------------------------------------------------------------------------    
     %%  Rainbow style colour maps
    
     case {'R1', 'R01', 'RAINBOW'}   % Reasonable rainbow colour map after it has been
                              % fixed by equalisecolourmap.
       desc = ['The least worst rainbow colour map I can devise.  Note there are' ...
              ' small perceptual blind spots at yellow and red'];
       attributeStr = 'rainbow';
       hueStr = 'bgyrm';       
       colpts = [35 63 -98
                 45 -14 -30
                 60 -55 60
                 85 0 80
                 55 60 62
                 75 55 -35];      
       sigma = 7;
       splineorder = 2;  % linear path 

     case {'R2', 'R02', 'RAINBOW2'}   % Similar to R1 but with the colour map finishing
                               % at red rather than continuing onto pink.
       desc = ['Reasonable rainbow colour map from blue to red.  Note there is' ...
              ' a small perceptual blind spot at yellow'];
       attributeStr = 'rainbow';
       hueStr = 'bgyr';
       colpts = [35 60 -98
                 45 -15 -30
                 60 -55 60
                 85 0 78
                 55 73 68];
       sigma = 5;
       splineorder = 2;  % linear path 
       
    
       
      case {'R2A', 'R02A'}   % 'Light', low contrast rainbow map for use with relief shading
        desc = 'Light, low contrast rainbow map for use with relief shading';
        attributeStr = 'rainbow';
        hueStr = 'bgyrm';
        colpts = [50  30 -80
                  60 -40 0
                  85 0 78
                  60 60 60
                  70 60 -20
                  80 40 -30];
        sigma = 7;
        splineorder = 2;  % linear path 
        W = [1,0,0];
        
        
      case {'R2AA', 'R02AA'}   % 'Light', low contrast rainbow map for use with relief shading
        desc = 'Light, low contrast rainbow map for use with relief shading';
        attributeStr = 'rainbow';
        hueStr = 'bgyrm';
        colpts = [50  30 -80
                  60 -40 0
                 85 0 78
                  60 60 60
                  70 60 -20];
        sigma = 7;
        splineorder = 2;  % linear path 
        W = [1,0,0];
        
        
      case {'R3', 'R03', 'RAINBOW3'}   % Diverging rainbow.  The blue and red
                                % points are matched in lightness and
                                % chroma as are the green and magenta
                                % points. Works well.
       desc = ['Diverging-rainbow colourmap. Yellow is the central reference' ...
               ' colour. The blue and red end points are matched in lightness' ...
               ' and chroma'];
       attributeStr = 'diverging-rainbow';
       hueStr = 'bgymr';
       colpts = [45 39 -83
                 52 -23 -23
                 60 -55 55
                 85 -2 85
                 60 74 -17
                 45 70 59];
       
       sigma = 5;
       splineorder = 2;  % linear path 
       
     % More vivid version of R2 with greater range of luminance. However I
     % find the vivid colours tend to make you segment regions in your data
     % on the basis of colour and you can lose the overall sense of structure.
     case {'R4', 'R04', 'RAINBOW4'}  
       desc = ['Rainbow colour map from blue to red.  Note there is' ...
              ' a small perceptual blind spot at yellow'];
       attributeStr = 'rainbow';
       hueStr = 'bgyr';
       colpts = [10 42 -57
                 20 59 -80
                 40 55 -93
                 50 -23 -25
                 60 -60 60
                 90 -12 89
                 70 33 72
                 55 75 65
                 45 70 55];
       sigma = 5;
       splineorder = 2;  % linear path 

       
       
       
       
     %-----------------------------------------------------------------------------    
     %%  Isoluminant colour maps       
      
     case {'I1', 'I01'}
      desc = ['Isoluminant blue to green to orange at lightness 70.  '...
              'Poor on its own but works well with relief shading'];
      attributeStr = 'isoluminant';
      hueStr = 'cgo';
      colpts = [70 ch2ab(40, -115)
                70 ch2ab(50, 160)
                70 ch2ab(50,  90)
                70 ch2ab(50,  45)];
      W = [1 1 1];

     case {'I2', 'I02'}  % Adaptation of I1 shifted to 80 from 70
      desc = ['Isoluminant blue to green to orange at lightness 80.  '...
              'Poor on its own but works well with relief shading'];
      attributeStr = 'isoluminant';
      hueStr = 'cgo';
      colpts = [80 ch2ab(36, -115)
                80 ch2ab(50, 160)
                80 ch2ab(50,  90)
                80 ch2ab(46,  55)];
      W = [1 1 1];
      
     case {'I3', 'I03'}
      desc = ['Isoluminant blue to pink at lightness 70.  '...
              'Poor on its own but works well with relief shading'];
      attributeStr = 'isoluminant';
      hueStr = 'cm';
      colpts = [70 ch2ab(40, -125)
                70 ch2ab(40, -80)
                70 ch2ab(40, -40)
                70 ch2ab(50,  0)];
      W = [1 1 1];      
      

     %-----------------------------------------------------------------------------    
     %% Colour Maps for the Colour Blind
     %% Protanopic / Deuteranopic colour maps
     
     % The Protanopic / Deuteranopic colour space is represented by two planes.
     % One plane is defined by the neutral axis and the colour point at yellow
     % 575nm, and the other defined by the neutral axis and the colour point at
     % blue 475nm.
     % Hue angles in a,b space of the key colours common to Protanopic/Deuteranopic
     % and Trichromatic viewers are:
     % yellow 575nm =  1.6165 radians,  92.62 degrees
     % blue   475nm = -1.3835 radians, -79.27 degrees

     case 'CBL1'  % Linear/diverging map for Protanopic/Deuteranopic viewers The
                  % symmetry requirements for diverging maps means that the
                  % colours are not as saturated as one would like. However it
                  % works better than CBL2
      desc = 'Linear map for Protanopic/Deuteranopic viewers';
      attributeStr = 'linear-protanopic-deuteranopic';
      hueStr = 'kbjyw';
      colpts = [5 0 0
                28 ch2ab(61, blue475)
                50 0 0
                72 ch2ab(61, yellow575) 
                95 0 0];
      W = [1 0 0];      
      splineorder = 3;
     
     case 'CBL2'  % Linear map with maximal chroma for Protanopic/Deuteranopic
                  % viewers. Does not work as well as one would hope.  The
                  % colours are too uneven.
      desc = 'Linear map with maximal chroma for Protanopic/Deuteranopic viewers';
      attributeStr = 'linear-protanopic-deuteranopic';
      hueStr = 'kbw';
      colpts = [5 0 0
                20 ch2ab(30, blue475)
                58 ch2ab(69, blue475)
                68 0 0
                85 ch2ab(98, yellow575)
                96 ch2ab(5, yellow575)                
                98 0 0];
      W = [1 0 0];    
      splineorder = 2;
      sigma = 9;       
      
     case 'CBL3'  % Linear map up the blue edge of the colour space
      desc = 'Linear blue map for Protanopic/Deuteranopic viewers';            
      attributeStr = 'linear-protanopic-deuteranopic';
      hueStr = 'kbw';
      colpts = [5 0 0
                56 ch2ab(68, blue475)
                95 0 0];
      W = [1 0 0];    
      splineorder = 2;
      sigma = 9;       
      
     case 'CBL4'  % Linear map up the yellow edge of the colour space
      desc = 'Linear yellow map for Protanopic/Deuteranopic viewers';
      attributeStr = 'linear-protanopic-deuteranopic';
      hueStr = 'kyw';
      colpts = [5 0 0
                25 ch2ab(34, yellow575)                
                88 ch2ab(88, yellow575)
                95 ch2ab(5, yellow575)];
      W = [1 0 0];      
      splineorder = 2;
      sigma = 9;      
      
     case 'CBD1'  % Diverging map blue-white-yellow
      desc = 'Diverging map blue-white-yellow';
      attributeStr = 'diverging-protanopic-deuteranopic';
      hueStr = 'bwy';
      colpts = [60 ch2ab(63, blue475)
                95 0 0
                60 ch2ab(63, yellow575)];
      W = [1 0 0];      
      splineorder = 2;
      sigma = 5;

     case 'CBD2'  % Diverging-linear map blue-grey-yellow
      desc = 'Diverging-linear map blue-grey-yellow';
      attributeStr = 'diverging-linear-protanopic-deuteranopic';
      hueStr = 'bjy';
      colpts = [57 ch2ab(67, blue475)
                73 0 0
                89 ch2ab(67, yellow575)];
      W = [1 0 0];      
      splineorder = 2;
      
     case 'CBC1'  % 4-phase cyclic map blue-white-yellow-black
      desc = '4-phase cyclic map blue-white-yellow-black';
      attributeStr = 'cyclic-protanopic-deuteranopic';
      hueStr = 'bwyk';
      colpts = [56 ch2ab(61, blue475)
                96 0 0
                56 ch2ab(61, yellow575)
                16 0 0];
      W = [1 0 0];      
      splineorder = 2;
      sigma = 5;
      
     case 'CBC2'  % 2-phase cyclic map white-yellow-white-blue
      desc = '2-phase cyclic map white-yellow-white-blue';
      attributeStr = 'cyclic-protanopic-deuteranopic';
      hueStr = 'wywb';
      colpts = [96 0 0
                55 ch2ab(68, yellow575)
                96 0 0
                55 ch2ab(68, blue475)];
      W = [1 0 0];      
      splineorder = 2;
      sigma = 5;
      
      
     %% Tritanopic colour maps
     % Hue angles in a,b space of the key colours common to Tritanopic
     % and Trichromatic viewers.
     % red  660nm  0.5648 radians,   32.36 degrees
     % cyan 485nm -2.4213 radians, -138.73 degrees
      
     case 'CBTL1'  % Tritanopic linear map with maximal chroma
      desc = 'Tritanopic linear map with maximal chroma';
      attributeStr = 'linear-tritanopic';
      hueStr = 'krjcw';
      colpts = [5 0 0
                20 ch2ab(50, red660)
                56 ch2ab(95, red660)
                78 ch2ab(52, cyan485)                
                98 0 0];
      W = [1 0 0];      
      splineorder = 2;
      sigma = 9;
     
     
     case 'CBTL2'  % Tritanopic linear map, can also be used as a diverging
                   % map. However the symmetry requirements of a diverging
                   % map results in colours of lower chroma
      desc = 'Tritanopic linear map';
      attributeStr = 'linear-tritanopic';
      hueStr = 'krjcw';
      colpts = [5 0 0
                25 ch2ab(58, red660)
                50 0 0
                75 ch2ab(58, cyan485)
                95 0 0];
      W = [1 0 0];      
      splineorder = 3;
      
     case 'CBTL3'  % Linear map up the blue green edge of the colour space
      desc = 'Tritanopic linear blue map';
      attributeStr = 'linear-tritanopic';
      hueStr = 'kcw';
      colpts = [5 0 0
                70 ch2ab(40, cyan485)
                85 ch2ab(40, cyan485)                
                95 0 0];
      W = [1 0 0];      
      splineorder = 3;

     case 'CBTL4'  % Linear map up the red edge of the colour space
      desc = 'Tritanopic linear red/heat map';
      attributeStr = 'linear-tritanopic';
      hueStr = 'krw';
      colpts = [5 0 0
                45 ch2ab(100, red660)
                70 ch2ab(50, red660)                
                95 0 0];
      W = [1 0 0];      
      splineorder = 3;
      
     case 'CBTD1'  % Tritanopic diverging map
      desc = 'Tritanopic diverging map';
      attributeStr = 'diverging-tritanopic';
      hueStr = 'cwr';
      colpts = [75 ch2ab(39, cyan485)
                98 0 0                
                75 ch2ab(39, red660)];
      W = [1 0 0];      
      splineorder = 2;
      sigma = 5;

     case 'CBTC1'  % 4-phase tritanopic cyclic map
      desc = '4-phase tritanopic cyclic map';
      attributeStr = 'cyclic-tritanopic';
      hueStr = 'cwrk';
      colpts = [70 ch2ab(39, cyan485)
                100 0 0                
                70 ch2ab(39, red660)
                40 0 0];
      W = [1 0 0];      
      splineorder = 2;
      sigma = 5;      

     case 'CBTC2'  % 2-phase tritanopic cyclic map
      desc = '2-phase tritanopic cyclic map';
      attributeStr = 'cyclic-tritanopic';
      hueStr = 'wrwc';
      colpts = [100 0 0                
                70 ch2ab(41, red660)
                100 0 0
                70 ch2ab(41, cyan485)];
      W = [1 0 0];      
      splineorder = 2;
      sigma = 5;      

      
     %----------------------------------------------------------------------------- 
     %%  Experimental colour maps and colour maps that illustrate some design principles


     case 'X1'
      desc = ['Two linear segments with different slope to illustrate importance' ...
              ' of lightness gradient.'];
      attributeStr = 'linear-lightnessnormalised';
      hueStr = 'by';
      colpts = [30 ch2ab(102, -54)
                40  0   0
                90  ch2ab(90, 95)];
      W = [1 0 0];
      splineorder = 2;
      
     case 'X2'
      desc = ['Two linear segments with different slope to illustrate importance' ...
              ' of lightness gradient.'];
      attributeStr = 'linear-CIE76normalised';
      hueStr = 'by';
      colpts = [30 ch2ab(102, -54)
                40  0   0
                90  ch2ab(90, 95)];
      W = [1 1 1];
      splineorder = 2;      

      
     case 'X3' % Constant lightness 'v' path to test unimportance of having a smooth
               % path in hue/chroma.  Slight 'feature' at the red corner (Seems more
               % important on poor monitors)
      attributeStr = 'isoluminant-HueChromaSlopeDiscontinuity';
      hueStr = 'brg';
      colpts = [50 17 -78
                50 77 57
                50 -48 50];
      splineorder = 2;  % linear path      
      W = [1 1 1];      
      
      
     % Yukky looking attempt at getting a constantly increasing lightness path
     case {'X5', 'X05'}
       desc = '';
       attributeStr = 'rainbow';
       hueStr = '';

       b20 = [20 59 -80];
       g40 = [40 -45 45];
       b40 = [40  25 -65];
       g50 = [50 -54 52];
       r55 = [55 80 66];
       o65 = [65 51 73];
       g80 = [80 -80 77];
       y90 = [90 -7 90];
       y97 = [97 -22 94];
       
       colpts = [b20
                 b40
                 r55
                 g80
                 y97];

       sigma = 5;
       splineorder = 2;  % linear path 

       
        % Rainbow constructed from an unwinding of the C6 cyclic
        % map. Interesting but I don't think it is the best thing to do
      case {'X6', 'X06'}  
        desc = 'Rainbow constructed from an unwinding of the C6 cyclic map'
        attributeStr = 'linear';
        hueStr = 'bcgyrm';
         
         y90 = [90 -7 90];
         m90 = [90 24 -17];
         c90 = [90 -48 -14];
         
         y85 = [85 -22 84];
         m85 = [85 35 -23];
         c85 = [85 -46 -14];
         
         y75 = [75 -16 74];
         m75 = [75 60 -40];
         c75 = [75 -39 -16];         
         
         r55 = [55 80 66];
         b55 = [55 19 -71];
         g55 = [55 -60 56];
         
         r50 = [50 78 62];
         b50 = [50 30 -81];
         g50 = [50 -54 52];
         
         colpts = [b50
                   c90
                   g50                   
                   y90
                   r50
                   m90];
              
         
         sigma = 7;
         splineorder = 2;  
         
       
       
       
       
       
       % A map inspired by Parula but with a smoother path that does not attempt to
       % include cyan and maintains a more uniform slope upwards in LAB space.
       % It works very much better than Parula!
      case 'X7'
        desc = 'Blue-Green-Orange-Yellow map similar to Parula but better!';
        attributeStr = 'linear';
        hueStr = 'bgoy';

        colpts = [25 50 -70
                  40 35 -75
                  55 -47 0 
                  70 -20 70
                  80  20 80
                  95 -21 92];
        sigma = 0;
        splineorder = 3;
       
       
       
     % A set of isoluminant colour maps only varying in saturation to test
     % the importance of saturation (not much) Colour Maps are linear with
     % a reversal to test importance of continuity.  
       
     case 'X10'  % Isoluminant 50 only varying in saturation
      attributeStr = 'isoluminant';
      hueStr = 'r';
      colpts = [50 0 0
                50 77 64
                50 0 0];
      splineorder = 2;
      sigma = 0;
      W = [1 1 1];
      
     case 'X11'  % Isoluminant 50 only varying in saturation
      attributeStr = 'isoluminant';
      hueStr = 'b';
      colpts = [50 0 0
                50 0 -56
                50 0 0];
      splineorder = 2;
      sigma = 0;
      W = [1 1 1];      
      
     case 'X12'  % Isoluminant 90 only varying in saturation
      attributeStr = 'isoluminant';
      hueStr = 'isoluminant_90_g';
      colpts = [90 0 0
                90 -76 80
                90 0 0];
      splineorder = 2;
      sigma = 0;
      W = [1 1 1];      

      
    % Difference in CIE76 and CIEDE2000 in chroma      
      
    case 'X13'  % Isoluminant 55 only varying in chroma. CIEDE76
      attributeStr= 'isoluminant-CIE76';
      hueStr = 'jr';
      colpts = [55 0 0
                55 80 67];
      splineorder = 2;
      W = [1 1 1];      
      formula = 'CIE76';
      
     case 'X14'  % Same as X13 but using CIEDE2000
      attributeStr= 'isoluminant-CIEDE2000';
      hueStr = 'jr';
      colpts = [55 0 0
                55 80 67];
      splineorder = 2;
      W = [1 1 1];      
      formula = 'CIEDE2000';      
      
     case 'X15'  % Grey 0 - 100. Same as No 1 but with CIEDE2000
      desc = 'Grey scale'; 
      attributeStr= 'linear-CIEDE2000';
      hueStr = 'grey';
      colpts = [  0 0 0      
                  100 0 0];
      splineorder = 2;
      formula = 'CIEDE2000';   

     case 'X16'  % Isoluminant 30 only varying in chroma
      attributeStr= 'isoluminant';
      hueStr = 'b';
      colpts = [30 0 0
                30 77 -106];
      splineorder = 2;
      W = [1 1 1];      
    
      
     case 'X21' % Blue to yellow section of rainbow map R1 for illustrating
                % colour ordering issues
       attributeStr= 'rainbow-section1';
       hueStr = 'bgy';       
       colpts = [35 60 -100
                45 -15 -30
                60 -55 60
                85 0 80];
       splineorder = 2;  % linear path 

     case 'X22' % Red to yellow section of rainbow map R1 for illustrating
                % colour ordering issues
       attributeStr= 'rainbow-section2';
       hueStr = 'ry';       
       colpts = [55 70 65
                85 0 80];
       splineorder = 2;  % linear path 

     case 'X23' % Red to pink section of rainbow map R1 for illustrating
                % colour ordering issues
       attributeStr= 'rainbow-section3';                
       hueStr = 'rm';       
       colpts = [55 70 65
                75 55 -35];      
       splineorder = 2;  % linear path 

       
     case 'XD1'  % Same as D1 but with no smoothing
      desc = 'Diverging blue-white-red colour map';
      attributeStr = 'diverging';
      hueStr = 'bwr';
      colpts = [40  ch2ab(83,-64)
                95  0   0
                40  ch2ab(83, 39)];   
      sigma = 0;
      splineorder = 2; 
      

     case 'X30'
      desc = 'red - green - blue interpolated in rgb';
      attributeStr = 'linear';
      hueStr = 'rgb';
      colourspace = 'RGB';
      colpts = [1.00 0.00 0.00
                0.00 1.00 0.00
                0.00 0.00 1.00];
      sigma = 0;
      W = [0 0 0];
      splineorder = 2; 
      
     case 'X31'
      desc = 'red - green - blue interpolated in CIELAB';
      attributeStr = 'linear';
      hueStr = 'rgb';
      colourspace = 'LAB';
      colpts = [53  80   67
                88 -86   83
                32  79 -108];
      sigma = 0;
      W = [0 0 0];
      splineorder = 2;       

      
     case 'XI4' % Maxiumum chroma isoluminant map at lightness 65
      colpts = [65 ch2ab(90, 135)
                65 ch2ab(88, 53)
                65 ch2ab(60, 10)
                65 ch2ab(90, -32)];
      W = [1 1 1];      
      splineorder = 2;
      sigma = 9;

     case 'XI5' % Maxiumum chroma isoluminant map at lightness 65
      colpts = [65 ch2ab(90, 135)
                65 ch2ab(90, -32)];
      W = [1 1 1];      
      splineorder = 2;

     case 'XI6' % Maxiumum chroma isoluminant map for deuteranope at
                % lightness 60
      yellow575 =  92.62; % degrees
      blue475   = -79.27; % degrees                
      colpts = [60 ch2ab(64, blue475)
                60 ch2ab(64, yellow575)];
      W = [1 1 1];      
      splineorder = 2;
      
      
      
      
     case 'XD7A'  % Linear diverging blue - magenta- grey - orange - yellow.
                 % Modified from 'D7' to have a double arch shaped path in an attempt
                 % to improve its Metric properties.  Also starts at lightness
                 % of 40 rather than 30.  The centre grey region is a bit too
                 % prominant and overall the map is perhaps a bit too 'bright'
      attributeStr = 'diverging-linear';
      hueStr = 'bmjoy';
      colpts = [40 ch2ab(88, -64)
                55 ch2ab(70, -30)
                64 ch2ab(2.5, -72.5)
                65 0 0       
                66 ch2ab(2.5, 107.5)
                75 ch2ab(70, 70)
                90 ch2ab(88,100)];
      splineorder = 3;
      
     case 'XD7C'  % Linear diverging  green - grey - yellow
                  % Reasonable, perhaps easier on the eye than D7

      attributeStr = 'diverging-linear';
      hueStr = 'gjy';
      rad = 65;
      colpts = [40 ch2ab(rad, 136)
                65 0 0                 
                90 ch2ab(rad, 95)];
      splineorder = 2;
      
      case 'XL102'  % Linear decreasing lightness, increasing saturation to blue
                    % Try to add intermediate colours YUK!
        attributeStr = 'linear';
        huestr = 'wb';
        colpts = [90  0  0
                  58  ch2ab(57, 130)
                  40  ch2ab(85, 45)
                  25  ch2ab(114, -54)];
        splineorder = 2;

      case 'XL103'  % Linear increasing lightness, increasing saturation to yellow
                    % Perhpaps the least objectionable increasing lightness 
                    % and saturation colour combination

        attributeStr = 'linear';
        huestr = 'km';
        colpts = [10  0  0
                  97 -22 94];
        splineorder = 2;

      case 'XL104'  % Linear increasing lightness, increasing saturation to magenta

        attributeStr = 'linear';
        huestr = 'km';
        colpts = [10  0  0
                  70 74 -46];
        splineorder = 2;
      
      
      
     %%-------------------------------------------------------------  
       
     otherwise
      
      % Invoke the catalogue search to help the user
      catalogue(I);
       
      clear map;
      clear name;
      clear desc;
      return
      
    end

    
    % Adjust chroma/saturation but only if colourspace is LAB
    if strcmpi(colourspace, 'LAB')
        colpts(:,2:3) = chromaK * colpts(:,2:3);
    end
    
    % Colour map path is formed via a b-spline in the specified colour space
    Npts = size(colpts,1);
    
    if Npts < 2
        error('Number of input points must be 2 or more')
    elseif Npts < splineorder
        splineorder = Npts;
        fprintf('Warning: Spline order is greater than number of data points\n')
        fprintf('Reducing order of spline to %d\n', Npts)
    end
    
    % Rely on the attribute string to identify if colour map is cyclic.  We may
    % want to construct a colour map that has identical endpoints but do not
    % necessarily want continuity in the slope of the colour map path.
    if strfind(attributeStr, 'cyclic')
        cyclic = 1;
        labspline = pbspline(colpts', splineorder, N);
    else
        cyclic = 0;
        labspline = bbspline(colpts', splineorder, N);
    end    
    
    % Apply contrast equalisation with required parameters. Note that sigma is
    % normalised with respect to a colour map of length 256 so that if a short
    % colour map is specified the smoothing that is applied is adjusted to suit.
    sigma = sigma*N/256;
    map = equalisecolourmap(colourspace, labspline', formula,...
                            W, sigma, cyclic, diagnostics);    

    % If specified apply a cyclic shift to the colour map
    if shift
        if isempty(strfind(attributeStr, 'cyclic'))
            fprintf('Warning: Colour map shifting being applied to a non-cyclic map\n');
        end
        map = circshift(map, round(N*shift));    
    end
    
    if reverse
       map = flipud(map);
    end    
        
    % Compute mean chroma of colour map
    lab = rgb2lab(map);
    meanchroma = sum(sqrt(sum(lab(:,2:3).^2, 2)))/N;
    
    % Construct lightness range description
    if strcmpi(colourspace, 'LAB')  % Use the control points
        L = colpts(:,1);
    else  % For RGB use the converted CIELAB values
        L = round(lab(:,1));
    end
    minL = min(L);
    maxL = max(L);
    
    if minL == maxL     % Isoluminant
        LStr = sprintf('%d', minL);
    
    elseif L(1) == maxL && ...
                 (~isempty(strfind(attributeStr, 'diverging')) ||...
                  ~isempty(strfind(attributeStr, 'linear')))
        LStr = sprintf('%d-%d', maxL, minL);
        
    else          
        LStr = sprintf('%d-%d', minL, maxL);
    end
    
    % Build overall colour map name
    name = sprintf('%s_%s_%s_c%d_n%d',...
                   attributeStr, hueStr, LStr, round(meanchroma), N);
    
    if shift
       name = sprintf('%s_s%d', name, round(shift*100)); 
    end
    
    if reverse
       name = sprintf('%s_r', name);
    end    
    
    
    if diagnostics  % Print description and plot path in colourspace
        fprintf('%s\n',desc);
        colourmappath(map, 'fig', 10)
    end
    
    
%------------------------------------------------------------------
% Conversion from (chroma, hue angle) description to (a*, b*) coords

function ab = ch2ab(chroma, angle_degrees)
    
    theta = angle_degrees/180*pi;
    ab = chroma*[cos(theta) sin(theta)];
    

%------------------------------------------------------------------
%
% Function to list colour maps with names containing a specific string.
% Typically this is used to search for colour maps having a specified attribute:
% 'linear', 'diverging', 'rainbow', 'cyclic', 'isoluminant' or 'all'.
% This code is awful!

function catalogue(str)
    
    if ~exist('str', 'var')
        str = 'all';
    end
    
    % Get all case expressions in this function
    caseexpr = findcaseexpr;
    
    % Construct all colour map names
    for n = 1:length(caseexpr)
        % Get 1st comma separated element in caseexpr
        label = strtok(caseexpr{n},',');  
        [~, name{n}] = cmap(label);
    end

    % Check each colour name for the specified search string.  Exclude the
    % experimental maps with label starting with X.
    fprintf('\n  CMAP label(s)            Colour Map name\n')
    fprintf('------------------------------------------------\n')
    
    found = 0;
    
    for n = 1:length(caseexpr)    
        if caseexpr{n}(1) ~= 'X'            
            if any(strfind(upper(name{n}), str)) || strcmpi(str, 'all')
                fprintf('%-20s   %s\n', caseexpr{n}, name{n});
                found = 1;
            end
        end
    end        

    if ~found
        fprintf('Sorry, no colour map with label or attribute %s found\n', str);
    end


%-----------------------------------------------------------------------
% Function to find all the case expressions in this function.  Yuk there must
% be a better way!  Assumes the case statement do not span more than one line

function caseexpr = findcaseexpr
    [fid, msg] = fopen([mfilename('fullpath') '.m'], 'r');
    error(msg);
    
    caseexpr = {};
    n = 0;
    line = fgetl(fid);

    while ischar(line) 
        [tok, remain] = strtok(line);
        
        if strcmpi(tok, 'case')
            % If we are here remain should be the case expression.
            % Remove any trailing comment from case expression
            [tok, remain] = strtok(remain,'%');            

            % Remove any curly brackets from the case expression
            tok(tok=='{') = [];
            tok(tok=='}') = [];
            tok(tok=='''') = [];

            n = n+1;
            caseexpr{n} = tok;            
        end
        line = fgetl(fid);
    end
    
    caseexpr = strtrim(caseexpr);
    
    fclose(fid);

%-----------------------------------------------------------------------
% Function to parse the input arguments and set defaults

function [I, N, chromaK, shift, reverse, diagnostics] = parseinputs(varargin)
    
    p = inputParser;

    numericORchar    = @(x) isnumeric(x) || ischar(x);
    numericORlogical = @(x) isnumeric(x) || islogical(x);
    
    % The first argument is either a colour map label string or a string to
    % search for in a colourmap name. If no argument is supplied it is assumed
    % the user wants to list all possible colourmaps.
    addOptional(p, 'I', 'all', @ischar); 
    
    % Optional parameter-value pairs and their defaults    
    % ** Note if you are using R2012 or earlier you should change the calls to
    % 'addParameter' below to 'addParamValue' and hopefully the code will work
    % for you.
    addParameter(p, 'N',     256, @isnumeric);  
    addParameter(p, 'shift',   0, @isnumeric);  
    addParameter(p, 'chromaK', 1, @isnumeric);     
    addParameter(p, 'reverse', 0, numericORlogical);  
    addParameter(p, 'diagnostics', 0, numericORlogical);  
    
    parse(p, varargin{:});
    
    I = strtrim(upper(p.Results.I));
    N = p.Results.N;
    chromaK     = p.Results.chromaK;
    shift       = p.Results.shift;
    reverse     = p.Results.reverse;
    diagnostics = p.Results.diagnostics;    
    
    if abs(shift) > 1
        error('Cyclic shift fraction magnitude cannot be larger than 1');
    end
    
    if chromaK < 0
        error('chromaK must be greater than 0')
    end    
    
    if chromaK > 1
        fprintf('Warning: chromaK is greater than 1. Gamut clipping may occur')
    end        
    
    
    
