function h=button(r,txt,cb);

% BUTTON:  Creates a GUI button.
%
% h=BUTTON(R,TXT,CB) defines a button with the position R, text TXT,
% and callback CB.  It uses UICONTROL function.
%

% svn $Id: button.m 436 2010-01-02 17:07:55Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                            A. Shcherbina          %
%===========================================================================%


h=uicontrol('units','normalized', ...
            'position',r, ...
	    'string',txt, ...
	    'callback',cb, ...
	    'interruptible','on');

if (nargout==0),
  clear;
end,

return


