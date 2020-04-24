function h=radiobox(r,tit,texts,name)

% RADIOBOX Creates a group of radiobuttons.
%
% h=RADIOBOX(r,tit,texts,name)
%  
% Creates a frame at position R, title TIT, and places a set of
% radiobuttons in it (one per cell in TEXTS). Some functionality is
% also implemented (e.g. exclusive selection). The global variable
% GUI will have a field NAME which indicates currently selected
% button, and a field NAME_h with the handles.
%

% svn $Id: radiobox.m 996 2020-01-10 04:28:56Z arango $
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

a=uicontrol('Style', 'frame',                                           ...
            'Units', 'normalized',                                      ...
            'Position', r);
r=r+[xspace yspace -2*xspace -2*yspace];
r(2)=r(2)+r(4)-dy;
r(4)=dy;

if (~isempty(tit))
  uicontrol('Style', 'text',                                            ...
            'String', tit,                                              ...
            'Units', 'normalized',                                      ...
            'Position', r,                                              ...
            'FontWeight', 'bold');
  r=next_line(r);
end

for k=1:L
  h(k)=uicontrol('Style', 'radiobutton',                                ...
                 'Units', 'normalized',                                 ...
                 'Position', r,                                         ...
                 'String', texts{k},                                    ...
                 'Value', k==1,                                         ...
                 'Callback',                                            ...
                 ['global GUI;'                                         ...
                  'GUI.' name '=' num2str(k) ';'                        ...
                  'for h=GUI.' name '_h;set(h,''value'',0); end;'       ...
                  'set(GUI.' name '_h(' num2str(k) '),''value'',1);']);
  r=next_line(r);
end

eval(['GUI.' name '_h=h;']);
eval(['GUI.' name '=1;']);

if (nargout==0)
  clear h
end

return


function r=next_line(r)

r(2)=r(2)-r(4);

return


