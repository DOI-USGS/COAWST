function h=textbox(r, tit, texts, name)

% TEXTBOX:  Creates a textbox in a frame
%
% h=TEXTBOX(r,tit,texts,name)
%  
% Creates a frame at position R, title TIT, and places a set of boxes
% in it (one per cell in TEXTS). The global variable GUI will have a
% field NAME_h with the handles.
%

% svn $Id: textbox.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                            A. Shcherbina        %
%=========================================================================%

global GUI

xspace=.01;
yspace=.01;
dy=.04;
L=length(texts);

a=uicontrol('style','frame',                                            ...
           'units','normalized',                                        ...
           'position',r);
r=r+[xspace yspace -2*xspace -2*yspace];
r(2)=r(2)+r(4)-dy;
r(4)=dy;

if (~isempty(tit))
  uicontrol('style','text',                                             ...
            'string',tit,                                               ...
            'units','normalized',                                       ...
            'position',r,                                               ...
            'fontweight','bold');
  r=next_line(r);
end

h(1:L)=0;
for k=1:L
  h(k)=uicontrol('style','edit',                                        ...
                 'units','normalized',                                  ...
                 'position',r,                                          ...
                 'string',texts{k});
  r=next_line(r);
end

eval(['GUI.' name '_h=h;']);

if (nargout==0)
  clear h;
end

return

function r=next_line(r)

r(2)=r(2)-r(4);

return


