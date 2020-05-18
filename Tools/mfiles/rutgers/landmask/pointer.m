function pointer(shape)

% POINTER:  Changes the mouse pointer.
%
% POINTER(shape) 
%
% Changes the pointer according to SHAPE.  If called without arguments,
% restores to normal arrow
%
% Available shapes:
%
%    AIM
%    POINT
%    RECT
%    ZOOM
%

% svn $Id: pointer.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                            A. Shcherbina        %
%=========================================================================%

% Set default arrow pointer.

if (nargin==0)
  set(gcf,'pointer','arrow');
  return
end

% Set various custum pointers.

switch lower(shape)
  case 'aim' 
    hs=[8,8];
    ptr=['ooooo++*++oooooo'
         'ooo++*****++oooo'
         'oo+**o+*+o**+ooo'
         'o+*ooo+*+ooo*+oo'
         'o+*oooo*oooo*+oo'
         '+*ooooo*ooooo*+o'
         '+*++ooooooo++*+o'
         '******ooo******o'
         '+*++ooooooo++*+o'
         '+*ooooo*ooooo*+o'
         'o+*oooo*oooo*+oo'
         'o+*ooo+*+ooo*+oo'
         'oo+**o+*+o**+ooo'
         'ooo++*****++oooo'
         'ooooo++*++oooooo'
         'oooooooooooooooo'];
  case 'zoom'
    hs=[7,7];
    ptr=['oooo****oooooooo'
         'oo**++++**oooooo'
         'o*++++++++*ooooo'
         'o*+++**+++*ooooo'
         '*++++**++++*oooo'
         '*++******++*oooo'
         '*++******++*oooo'
         '*++++**++++*oooo'
         'o*+++**+++*ooooo'
         'o*++++++++**oooo'
         'oo**++++*****ooo'
         'oooo****ooo***oo'
         'oooooooooooo***o'
         'ooooooooooooo***'
         'oooooooooooooo*o'
         'oooooooooooooooo'];
  case 'rect'
    hs=[5,5];
    ptr=['oooo*ooooooooooo'
         'oooo*+oooooooooo'
         'oooo*+oooooooooo'
         'oooo*+oooooooooo'
         '****o****ooooooo'
         'o+++*+++++oooooo'
         'oooo*+oooooooooo'
         'oooo*+oooooooooo'
         'oooo*+oo******oo'
         'ooooo+oo*++++*oo'
         'oooooooo*++++*oo'
         'oooooooo******oo'
         'oooooooooooooooo'
         'oooooooooooooooo'
         'oooooooooooooooo'
         'oooooooooooooooo'];
  case 'point'
    hs=[9,3];
    ptr=['oooooooooo*+oooo'
         'ooooooooo*+ooooo'
         'oooooooo*+oooooo'
         'oo*+ooo*+ooooooo'
         'oo*+oo*+oooooooo'
         'oo*+o*+ooooooooo'
         'oo*+*+oooooooooo'
         'oo**+ooooooooooo'
         'oo*+oooooooooooo'
         'oooooooooooooooo'
         'oooooooooooooooo'
         'oooooooooooooooo'
         'oooooooooooooooo'
         'oooooooooooooooo'
         'oooooooooooooooo'
         'oooooooooooooooo'];
  otherwise
    set(gcf,'pointer',shape);
    return
end

ptr=double(ptr);
ptr(ptr=='*')=1;
ptr(ptr=='+')=2;
ptr(ptr=='o')=NaN;
set(gcf,'pointer','custom',                                             ...
        'pointershapecdata',ptr,                                        ...
        'PointerShapeHotSpot',hs)

return
