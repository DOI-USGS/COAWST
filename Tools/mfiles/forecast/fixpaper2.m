%   FIXPAPER  Makes the hard copy look like the screen.
%   Sets the 'paperposition' so that what goes onto paper 
%   (or into the postscript file), looks just like what 
%   is on your computer screen.

%   All we do here is set the width and height of the 
%   'paperposition' to be the same as the width and height
%   of the current figure on the screen.

%   Rich Signell
un=get(gcf,'units');
pun=get(gcf,'paperunits');
set(gcf,'units','inches');
pos=get(gcf,'pos');

% if figure window window > height, use landscape
%if(pos(3)>pos(4)),
%
%else 
%  orient('portrait')
%end

% set the paperposition to the same as the figure window,
% starting 1/2 inch from the bottom left hand corner of the page.

set(gcf, 'paperunits', 'inches', 'paperposition',[.5 .5 pos(3:4)]);

%reset paper units to original values
set(gcf,'units',un);
set(gcf,'paperunits',pun);
