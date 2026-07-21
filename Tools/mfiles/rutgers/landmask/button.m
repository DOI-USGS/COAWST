function h=button(r,txt,cb);

% BUTTON:  Creates a GUI button.
%
% h=BUTTON(R,TXT,CB)
%  
% defines a button with the position R, text TXT, and callback CB.
% It uses UICONTROL function.
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                             A. Shcherbina        %
%=========================================================================%


h=uicontrol('units','normalized',                                       ...
            'position',r,                                               ...
	    'string',txt,                                               ...
	    'callback',cb,                                              ...
	    'interruptible','on');

if (nargout==0)
  clear;
end

return


