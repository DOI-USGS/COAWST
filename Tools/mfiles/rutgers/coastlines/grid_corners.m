%
%  GRID_CORNERS:  Driver to set the grid corners for "seagrid"
%
%  This is a user modifiable script that can be used to set the grid corners
%  for "seagrid" to guarantee rectangular grids.  Such values can written
%  into a boundary.dat which can be read by "seagrid".
%

% svn $Id: grid_corners.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           John Wilkin             %
%===========================================================================%

location = 'eauc';
location = 'nena'

switch location
  
  case 'eauc' % east Auckland current

    % axis limits
    ax = [170 182 -40 -30];

    % projection
    m_proj('mercator','lon',ax(1:2),'lat',ax(3:4));

    [xakl,yakl]     = m_ll2xy(174.84,-36.87);          % Auckland 
    [xecape,yecape] = m_ll2xy(178.53,-37.63);          % East Cape
    [xn1,yn1]       = m_ll2xy(171,-33.5,'clip','off'); % northeast 1000 m dep
    [xker,yker]     = m_ll2xy(172,-34.1);              % Kermadec Is.
    [xse1,yse1]     = m_ll2xy(179.4,-37);              % southeast limit
    [xn2,yn2]       = m_ll2xy(171,-33.5,'clip','off'); % 

    % calculate the corners of the desired rectangular box
    
    % specify nw corner
    [xnw,ynw] = m_ll2xy(171.7,-33.8,'clip','off');
    
    % fit a line through nw point and akl: y=m*x+c
    A = [xnw 1; xakl 1];
    b = [ynw; yakl];
    tmp = A\b;
    m1 = tmp(1);
    c1 = tmp(2);

    % get coords for sw corner 
    xsw = m_ll2xy(177.3,0,'clip','off'); % lat doesn't matter in mercator
    ysw = m1*xsw + c1;
    
    % fit a line through nw point perpendicular to previous line
    m2 = -1/m1;
    c2 = ynw - m2*xnw;
    
    % get coords for ne corner 
    xne = m_ll2xy(175.1,0,'clip','off'); % lat doesn't matter in mercator
    yne = m2*xne + c2;
    
    % get se corner
    xse = xne + (xsw-xnw);
    yse = yne + (ysw-ynw);

  case 'nena' % Northeast North America

    % axis limits
    ax = [265 315 16 55];

    % projection
    m_proj('mercator','lon',ax(1:2),'lat',ax(3:4));

    % set coastal locations to define western boundary 
    [x2,y2]   = m_ll2xy(285.9,40.56,'clip','off');  % Staten island
    [x2,y2]   = m_ll2xy(283.65,39.18,'clip','off'); % upper Chesapeake Bay
    [x2,y2]   = m_ll2xy(289.9,43.7,'clip','off'); % Gulf of Maine coast
    [x1,y1]   = m_ll2xy(273.6,30.46,'clip','off'); % Niceville
    [x1,y1]   = m_ll2xy(272.0,30.69,'clip','off'); % Mobile
    
    % fit a line through the points corners
    A = [x1 1; x2 1];
    b = [y1; y2];
    tmp = A\b;
    m1 = tmp(1);
    c1 = tmp(2);

    % get coords for nw corner
    nw_long = 300.0; % longitude 50W
    xnw = m_ll2xy(nw_long,0,'clip','off'); % lat doesn't matter in mercator
    ynw = m1*xnw + c1;

    % get coords for sw corner
    sw_long = 270.0; % longitude 90W
    sw_long = 269.6; % longitude approx 90W
    xsw = m_ll2xy(sw_long,0,'clip','off'); % lat doesn't matter in mercator
    ysw = m1*xsw + c1;
    
    % fit a line through nw point perpendicular to previous line
    m2 = -1/m1;
    c2 = ynw - m2*xnw;
    
    % get coords for ne corner 
    ne_long = 311.0;
    ne_long = 310.0;
    ne_long = 309.0;
    xne = m_ll2xy(ne_long,0,'clip','off'); % lat doesn't matter in mercator
    yne = m2*xne + c2;
    
    % get se corner
    xse = xne + (xsw-xnw);
    yse = yne + (ysw-ynw);

    if 0
    % set these two corners
    [xsw,ysw]   = m_ll2xy(272.47,30.0,'clip','off');  % sw corner mid GoMex
    [xnw,ynw]   = m_ll2xy(297.61,50.32,'clip','off'); % nw corner Canada

    % calculate the other corners of the desired rectangular box
    
    % fit a line through western corners
    A = [xnw 1; xsw 1];
    b = [ynw; ysw];
    tmp = A\b;
    m1 = tmp(1);
    c1 = tmp(2);

    % fit a line through nw point perpendicular to previous line
    m2 = -1/m1;
    c2 = ynw - m2*xnw;
    
    % get coords for ne corner 
    ne_long = 308.0;
    xne = m_ll2xy(ne_long,0,'clip','off'); % lat doesn't matter in mercator
    yne = m2*xne + c2;
    
    % get se corner
    xse = xne + (xsw-xnw);
    yse = yne + (ysw-ynw);
    end

end

% Convert to lon/lat, using ROMS corner numbering convention

[lon_1,lat_1] = m_xy2ll(xnw,ynw);
[lon_2,lat_2] = m_xy2ll(xsw,ysw);
[lon_3,lat_3] = m_xy2ll(xse,yse);
[lon_4,lat_4] = m_xy2ll(xne,yne);
corners.lon = [lon_1 lon_2 lon_3 lon_4]';
corners.lat = [lat_1 lat_2 lat_3 lat_4]';

% Quick plots

figure(1)
clf
plot(corners.lon([1:4 1]),corners.lat([1:4 1]),'-') 
amerc

figure(2)
clf
han = m_line(corners.lon([1:4 1]),corners.lat([1:4 1]));
set(han,'color','r')
m_grid
switch location
  case 'eauc'
    m_coast
  case 'nena'
    cst = load('natl_coast_int');
    han = m_line(cst.lon+360,cst.lat);
    set(han,'color','b');
end

% Create a corners files for seagrid

if min(corners.lon)>180
  corners.lon = corners.lon-360;
end
if max(corners.lon>180)
  warning([ 'Range of longitudes is ' mat2str(range(corners.lon))])
end
corners_data = [corners.lon corners.lat ones([4 1])];

% Save the boundary.dat file for seagrid

savefile = [ 'boundary_' location '.dat'];
if exist(savefile)==2
  reply = input([savefile ' exists. Overwrite? (y/n) '],'s');
  if strcmp(lower(reply),'y')
    save(savefile,'corners_data','-ascii')
    disp([ 'Wrote ' savefile])
  end
else
  save(savefile,'corners_data','-ascii')
  disp([ 'Wrote ' savefile])
end
