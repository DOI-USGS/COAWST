function h=axisscroll(type);

% AXISSCROLL Creates scrollbar(s) for the current axes.
%
% h=AXISSCROLL(TYPE) place a scrollbar according with the value of
% TYPE: 'l' left axis, 'r' right axis, 't' top axis, or 'b' bottom
% axis. If the axis limits are changed, like zooming in, a call to
% AXISSCROLL without arguments will allow to move the window. When
% zooming out to original limits the scrollbar(s) are hidden.
%
% Figure resizing is handeled automatically.
%
% Example:
%
%   pcolor(toeplitz(1:20))
%   axisscroll('t');
%   axisscroll('r');
%   xlim([3 6]);
%   ylim([3,6]);
%   axisscroll
%

% svn $Id: axisscroll.m 436 2010-01-02 17:07:55Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                            A. Shcherbina          %
%===========================================================================%

TAG='AxisScroll';
if (nargin==0),                          % just update positions
  xx=findobj(gcf,'tag',TAG);
  if (isempty(xx)),
    return;
  end,
  for hh=xx',
    positionscroll(hh);
  end,
else,
  UD.type=type;
  UD.xlim=xlim;
  UD.ylim=ylim;
  xx=findobj(gcf,'tag',TAG);             % delete xx
  h=uicontrol('style','slider', ...
              'tag',TAG, ...
	      'userdata',UD);
  positionscroll(h);
  set(gcf,'resizefcn','axisscroll');
  if (nargout==0)
    clear;
  end,
end

function positionscroll(h)

dist=2;
width=12;
ud=get(h,'userdata');
type=ud.type;

ou=get(gca,'units');
set(gca,'units','points');
r=get(gca,'position');
set(gca,'units',ou);

% Horizontal scroll, bottom of the figure.

if (type(1)=='b'),
  set(h,'units','points', ...
        'position',[r(1) r(2)-dist-width r(3) width], ...
        'callback', ['xl=xlim; p=get(gcbo,''value'');' ...
                     'xlim([0 diff(xl)]+p);']);

% Horizontal scroll, top of the figure.

elseif (type(1)=='t'),

   set(h,'units','points', ...
         'position',[r(1) r(2)+r(4)+dist r(3) width], ...
         'callback', ['xl=xlim; p=get(gcbo,''value'');' ...
                      'xlim([0 diff(xl)]+p);']);

% Vertical scroll bar, left of the figure.

elseif (type=='l'),

  set(h,'units','points', ...
        'position',[r(1)-dist-width r(2) width r(4)], ...
        'callback',['yl=ylim;p=get(gcbo,''value'');' ...
                    'ylim([0 diff(yl)]+p);']);

% Vertical scroll bar, right of the figure.

elseif (type=='r'),

  set(h,'units','points', ...
         'position',[r(1)+r(3)+dist r(2) width r(4)], ...
         'callback',['yl=ylim; p=get(gcbo,''value'');' ...
                     'ylim([0 diff(yl)]+p);'] );
end,

% Set limits.

if (type=='l' | type=='r'),
  yl=ylim;
  yl0=ud.ylim;
  set(h,'min',min(yl0)-1e-5, ...
        'max',max(yl0)-diff(yl), ...
        'value',min(yl), ...
        'sliderstep',[.1 1]*diff(yl)./diff(yl0));
   if (all(yl==yl0)),
     set(h,'visible','off');
   else,
     set(h,'visible','on');
   end,
else,
  xl=xlim;
  xl0=ud.xlim;
  set(h,'min',min(xl0)-1e-5, ...
        'max',max(xl0)-diff(xl), ...
	'value',min(xl), ...
	'sliderstep',[.1 1]*diff(xl)./diff(xl0));
  if (all(xl==xl0)),
    set(h,'visible','off');
  else,
    set(h,'visible','on');
  end,
end,

return
