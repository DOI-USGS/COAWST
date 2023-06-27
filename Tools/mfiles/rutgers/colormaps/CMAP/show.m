% SHOW - Displays an image with the right size, colors, range, and with a title.
%
% Usage:   
%         h = show(im);
%         h = show(im, figNo, title, clim, colourmap)
%
% Arguments:  im    - Either a 2 or 3D array of pixel values or the name
%                     of an image file;
%             figNo - Optional figure number to display image in. If
%                     figNo is 0 the current figure or subplot is
%                     assumed.  The default is to create a new figure.
%             title - Optional string specifying figure title. 
%                     Defaults to the variable name of the image.
%             clim  - Optional 2-vector specifying range of values to
%                     display. Useful for dealing with outlying values im
%                     your data, or for displaying data with a diverging
%                     colour map properly. Defaults to the full data range.
%         colourmap - Optional Nx3 matrix specifying a colour map.
%                     Defaults to gray(256).
%
% Returns:    h     - Handle to the figure.  This allows you to set
%                     additional figure attributes if desired.
%
% Apart from 'im' all arguments are optional and can be supplied in any
% order.  They are inferred from their types.
%
% Where possible the image is displayed as 'TrueSize', that is, pixels on the
% screen match pixels in the image.
%
% Unless you are doing a subplot (figNo==0) the window is sized to match
% the image, leaving no border,  hence saving desktop real estate.
%
% See also: COLORCET, SHOWSURF

% Copyright (c) 2000-2018 Peter Kovesi
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

% October   2000  Original version
% March     2003  Mods to alow figure name in window bar and allow for subplots.
% April     2007  Proper recording and restoring of MATLAB warning state.
% September 2008  Octave compatible
% May       2009  Reworked argument handling logic for extra flexibility
% January   2013  More Octave compatibility and proper restoring of warning state (Carnë Draug)
% March     2013  Ensure grey colour map has 256 levels rather than the
%                 default 64
% October   2017  Let imagesc apply default scaling for images < 500 pixels in size.
% April     2018  Improved argument handling.  Allow colour map and display
%                 range limits to be specified.

function h = show(im, varargin)

    Octave = exist('OCTAVE_VERSION', 'builtin') == 5; % Are we running under Octave

    s = warning('query','all');                 % Record existing warning state.
    warn_state = onCleanup (@() warning(s));    % Restore warnings at the end
    warning('off');                             % Turn off warnings that might arise if image
                                                % has to be rescaled to fit on screen

    [figNo, Title, clim, colourmap] = checkargs(varargin(:));
    
    % Check case where im is an image filename rather than image data
    if ~isnumeric(im) && ~islogical(im) 
        if isempty(Title)
            Title = im;        % Use file name for title
        end
        im = imread(im);
    elseif isempty(Title)
        Title = inputname(1);  % Use variable name of image data for title
    end

    sze = max(size(im));       % Maximum dimension of image

    if figNo > 0               % We have a valid figure number
        figure(figNo);         % Reuse or create a figure window with this number
        subplot('position',[0 0 1 1]); % Use the whole window
    elseif figNo == -1
        figNo = figure;        % Create new figure window
        subplot('position',[0 0 1 1]); % Use the whole window
    end

    if ndims(im) == 2          % Apply colour map
        if isempty(clim)
            imagesc(im);
        else
            imagesc(im, clim);
        end
        
        colormap(colourmap);  
    else
        imshow(im(:,:,1:3));   % Display as RGB (ignore any alpha channel)
    end

    if figNo == 0              % Assume we are trying to do a subplot 
        figNo = gcf;           % Get the current figure number
        axis('image'); axis('off');
        title(Title);          % Use a title rather than rename the figure
    else
        axis('image'); axis('off');
        set(figNo,'name', ['  ' Title])

        % If not running Octave and size of image > 500 plot image at 1:1
        % resolution. Otherwise we let imagesc use its default scaling.
        if ~Octave && sze > 500
            truesize(figNo);
        end
    end

    if nargout == 1
       h = figNo;
    end

    
%----------------------------------------------------------------------
%
% Process optional arguments. If an arguments is:
% - a scalar assume it is the figure number
% - a string, assume it is the figure title
% - a 1x2 vector assume it is the range display limits
% - a Nx3 matrix assume it is a colour map

function [figNo, Title, clim, colourmap] = checkargs(arg)
    
    % Set defaults
    figNo = -1;  % Default indicates create new figure
    Title = '';
    clim = [];
    colourmap = gray(256);

    % Overwrite defaults with anything we recognize in the arguments
    for n = 1:length(arg);
        if isscalar(arg{n})
            figNo = arg{n};
            
        elseif ischar(arg{n})
            Title = arg{n};
            
        elseif isnumeric(arg{n}) && all(size(arg{n}) == [1,2])
            clim = arg{n};
            
        elseif isnumeric(arg{n}) && size(arg{n},1) > 1 && size(arg{n},2) == 3
            colourmap = arg{n};            
            
        else
            error('Unable to process arguments');
        end
    end
    