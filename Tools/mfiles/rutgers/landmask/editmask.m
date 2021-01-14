function editmask(grid_file, varargin)

% EDITMASK:  Interactive Land/Sea mask editing for ROMS
%
% editmask(grid_file)
% editmask(grid_file, 'f')                    GSHHS
% editmask(grid_file, 'mycoast.mat')          see "get_coast.m"
% editmask(grid_file, 'ijcoast.mat')          see "ijcoast.m"
% editmask(grid_file, clon, clat)             see "get_coast.m"
%
% GUI for manual editing of ROMS Land/Sea mask on RHO-points.  To
% accelerate the processing, the Land/Sea mask is edited in (I,J)
% grid coordinates. If the (I,J) coordinates are not provided, it
% computes and writes them into file. If called without, one or
% both arguments, it will prompt for the needed file name(s). If
% the coastline data is in grid_file, it will read it and convert
% it to (I,J) fractional coordinates.
%
% Warning: 'editmask' is called recursive to process various Land/Sea
%          editing tasks.
%
% On Input:
%
%    grid_file   ROMS Grid NetCDF file name containing the grid
%                and mask arrays (string)
%
%            or, recursively editing tasks values: 'zoomin', 'zoomout',
%                'click', 'refresh', 'move', 'save', 'undo', 'exit' 
%
%    database    GSHHS database (character, OPTIONAL)
%                  'f'      full resolution
%                  'h'      high resolution (default)
%                  'i'      intermediate resolution
%                  'l'      load resolution
%                  'c'      crude resolution
%
%            or, a matlab (*.mat) containing coasline data
%                  (string, OPTIONAL)
%
%                  lon      coastline longitude
%                  lat      coastline latitude
%
%            or, C.Icst     coastline I-fractional indices
%                C.Jcst     coastline J-fractional indices
%                C.lon      coastline longitude
%                C.lat      coastline latitude
%
% or
%    clon        Coastline longitude (real vector, OPTIONAL)
%
%    clat        Coastline latitude  (real vector, OPTIONAL)
%
% Mouse shortcuts:
%
%    double click ==> Zoom in
%    right  click ==> Zoom out
%    middle click ==> change editing mode
%

%    Calls: WRITE_MASK and UVP_MASKS functions.
%           BUTTON, RADIOBOX, TEXTBOX, AXISSCROLL,
%

% svn $Id: editmask.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                            A. Shcherbina        %
%=========================================================================%

% Define and initialize persistent variables.

persistent changed rmask rlon rlat mask hplot Lp Mp mx my
persistent mfile fig zooming xx yy xl yl xcst ycst

% Single global variable to pass info to/from the callback routines.

global GUI

FIGURE_NAME='Land/Sea Mask Editor';

% Set colormap for land and sea.

 land1=[1.0000 0.9254 0.5451]; land2=[0.8039 0.7451 0.4392];  % goldenrod
%land1=[0.9333 0.8353 0.7176]; land2=[0.8039 0.7176 0.6196];  % bisque
%land1=[0.4745 0.3765 0.2980]; land2=[0.6706 0.5841 0.5176];  % brown

%seas1=[0.0000 0.9000 1.0000]; seas2=[0.0000 1.0000 1.0000];  % blue
%seas1=[0.0000 1.0000 1.0000]; seas2=[0.0000 0.8039 0.8039];  % cyan
 seas1=[0.2509 0.8784 0.8157]; seas2=[0.0000 0.7725 0.8039];  % turquoise

CMAP=[land1; land2; seas1; seas2];

% Set coastline line color and width.

LineColor='k';
LineWidth=1;

% Check input GRID file argument.

if (nargin < 1 || isempty(grid_file))
  [fn,pth]=uigetfile('*.nc','Select ROMS grid file...');
  if (~fn)
    return
  end
  grid_file=[pth,fn];
end
   
%==========================================================================
% The EDITMASK function is also used as a CallBack for some of the
% uicontrols,  so we need to figure out, what's required by looking
% at the 'grid_file' argument. Notice that EDITMASK is called recursively
% and 'grid_file' is the editing task.
%==========================================================================

switch lower(grid_file)

%--------------------------------------------------------------------------
% Zoom-in.
%--------------------------------------------------------------------------

  case 'zoomin'

    disable_click;
    editmask move;
    zooming=1;
    waitforbuttonpress
    xx0=xx; yy0=yy;                         % save current pointer position
    rbbox;                                  % track rubberband rectangle
    xx1=xx; yy1=yy;                         % new pointer position
    if ((xx0 ~= xx1) && (yy0 ~= yy1))       % trim limits and set zoom
      xx0(xx0<0)=0;  xx0(xx0>xl(2))=xl(2);
      xx1(xx1<0)=0;  xx1(xx1>xl(2))=xl(2);
      yy0(yy0<0)=0;  yy0(yy0>yl(2))=yl(2);
      yy1(yy1<0)=0;  yy1(yy1>yl(2))=yl(2);

      xlim([min([xx0 xx1]), max([xx0 xx1])]);
      ylim([min([yy0 yy1]), max([yy0 yy1])]);
      axisscroll;
    end
    enable_click;
    zooming=0;

%--------------------------------------------------------------------------
% Zoom-out.
%--------------------------------------------------------------------------

  case 'zoomout'

    xlim(xl);
    ylim(yl);
    axisscroll;

%--------------------------------------------------------------------------
% Edit Land/Sea mask.
%--------------------------------------------------------------------------

  case 'click'

    button=get(gcf, 'SelectionType');
    if (strcmp(button,'alt'))                 % zoom out on right click
      editmask zoomout;
      return
    end
    if (strcmp(button,'open'))                % zoom in on double click
      editmask zoomin;
      return
    end
    if (strcmp(button,'extend'))              % cycle modes on middle click
      m=mod(GUI.mode,3)+1;
      eval(get(GUI.mode_h(m),'callback'));
      editmask zoomout;
      return
    end
    if (within(xx,xlim) && within(yy,ylim))   % left click within edit area
      disable_click;
      switch GUI.tool
        case 1                                % point edit
          ix=floor(xx+1.5);
          iy=floor(yy+1.5);
          switch GUI.mode
            case 1                            % toggle between land and sea
              rmask(ix,iy)=~rmask(ix,iy);
            case 2                            % set land
              rmask(ix,iy)=0;
            case 3                            % set sea
              rmask(ix,iy)=1;
          end
        case 2                                  % area edit: save current
          xx0=xx; yy0=yy;                       % pointer position, track
          rbbox;                                % rubberband rectangle new
          xx1=xx; yy1=yy;                       % pointer position
          idsel=find((mx-xx0-1).*(mx-xx1-1)<=0 &                        ...
                     (my-yy0-1).*(my-yy1-1)<=0);% indices within rectangle
          switch GUI.mode
            case 1
              rmask(idsel)=~rmask(idsel);       % toggle between land/sea
            case 2
              rmask(idsel)=0;                   % set land
            case 3
              rmask(idsel)=1;                   % set sea
          end
      end
      changed=1;
      enable_click;
      editmask refresh;
    end

%--------------------------------------------------------------------------
% Update mask by changing colormap.
%--------------------------------------------------------------------------

  case 'refresh'

    cdata=rmask*2+mod(mod(mx,2)+mod(my,2),2)+1;
    set(hplot,'cdata',cdata');
    nm=[FIGURE_NAME,':  ', mfile];
    if (changed)
      nm=[nm,', modified...'];
    end
    set(fig,'name',nm);

%--------------------------------------------------------------------------
% Pointer movement: update pointer coordinates on menu.
%--------------------------------------------------------------------------

  case 'move'

    xy=get(gca,'currentpoint');
    xx=xy(1,1); yy=xy(1,2);
    if (within(xx,xlim) && within(yy,ylim))
      s=sprintf('[%d,%d]',floor(xx+.5),floor(yy+.5));
      if (zooming)
        pointer('zoom');
      elseif (GUI.tool==1)
        pointer('crosshair');
      elseif (GUI.tool==2)
        pointer('rect');
      end
    else
      s='---';
      pointer('arrow');
    end
    set(GUI.pos_h,'string',s);

%--------------------------------------------------------------------------
% Compute U-, V- and PSI masks.  Write out into GRID NetCDF file.
%--------------------------------------------------------------------------

  case 'save'

    [umask,vmask,pmask]=uvp_masks(rmask);
    mask=rmask;
    write_mask(mfile,rmask,umask,vmask,pmask);
    changed=0;
    editmask refresh;

%--------------------------------------------------------------------------
% Undo changes: restore last saved mask.
%--------------------------------------------------------------------------

  case 'undo'

    rmask=mask;
    changed=0;
    editmask refresh;

%--------------------------------------------------------------------------
% Done: inquire to save mask changes.
%--------------------------------------------------------------------------

  case 'exit'

    if (~changed)
      delete(gcf);
    else
      res=questdlg('The Land/Sea mask was changed, save?',FIGURE_NAME);
      switch res
        case 'Yes'
          editmask save;
          disp('Land/Sea mask was saved');
          delete(gcf);
        case 'No'
          disp('Land/Sea mask was NOT saved');
          delete(gcf);
      end
    end

%--------------------------------------------------------------------------
% Initialize: read mask and coastline data.
%--------------------------------------------------------------------------

  otherwise

% Kill all windows.

    delete(findobj(0,'tag','MaskEditor'));

% Process optional arguments.

    disp(' ')
    disp(['Processing grid file: ', grid_file]);
    tic;

    [pdir,fname,ext]=fileparts(grid_file);

    got_clon=false; got_clat=false;    % provided coastline data vectors
    extract_coast=false;               % extract coastlines from GSHHS
    got_matfile=false;                 % coastline from provided .mat file

    Cfile=which('gshhs_h.b','-ALL');   % select first directory found

    switch numel(varargin)
      case 0
        if (~isempty(Cfile))
          DIR=fileparts(Cfile{1});            % others are shadowed
          database='h';
          Cname=fullfile(DIR, 'gshhs_h.b');   % high resolution GSHHS
          extract_coast=true;                 % (default)
        else
          error('Cannot find GSHHS coastline dataset');
        end
      case 1
        if (ischar(varargin{1}))
          if (strfind(varargin{1}, '.mat'))
            coast_file=varargin{1};
            Cname=coast_file;
            if (~exist(coast_file,'file'))
              error(['Cannot file: ', coast_file]);
            end
            got_matfile=true;
          else
            database=varargin{1};
            switch database
              case 'f'                        % full resolution
                name='gshhs_f.b';
              case 'h'                        % high resolution
                name='gshhs_h.b';
              case 'i'                        % intermediate resolution
                name='gshhs_i.b';
              case 'l'                        % low resolution
                name='gshhs_l.b';
              case 'c'                        % crude resolution
                name='gshhs_c.b';
              otherwise
                error(['illegal GSHHS dataset resolution, ',database])
            end
%           if (~isempty(DIR))
            if (~isempty(Cfile))
              DIR=fileparts(Cfile{1});           % others are shadowed
              Cname=fullfile(DIR, name);      % selected GSHHS resolution
              extract_coast=true;
            else
              error(['Cannot find directory for: ', name]);
            end
          end
        elseif (isnumeric(varargin{1}))
          C.lon=varargin{1};                  % provide coastline longitude
          got_clon=true;
        end
      case 2
        C.lat=varargin{2};
        got_clat=true;
    end

% Set GRID structure.

    G=get_roms_grid(grid_file);
    mfile=grid_file;

% Get data from GRID structure.

    got_coast=all(isfield(G,{'lon_coast','lat_coast'}));  % GRID NetCDF has
                                                          % coastline data
    spherical=G.spherical;
    if (spherical)
      rlon=G.lon_rho;
      rlat=G.lat_rho;
    else
      rlon=G.x_rho;
      rlat=G.y_rho;
    end
    mask=G.mask_rho;
    rmask=mask;

    [Lp,Mp]=size(mask);
    [mx,my]=ndgrid(1:Lp,1:Mp);
    toc;

% Process coastline data.

    if (spherical)

      disp(' ');
      disp(['Processing coastline file: ', Cname]);
      tic;
    
      if (extract_coast)            % extract coastline from GSHHS dataset

        ijcoast=strcat(fname, '_ijcst_', database, '.mat');

        if (~exist(ijcoast, 'file'))
          dx=5*abs(mean(mean(diff(rlon))));
          dy=5*abs(mean(mean(diff(rlat))));
          if (dx == 0)
            dx=1.5;
          end
	  if (dy == 0)
            dy=1.5;
          end

          Llon=min(rlon(:));  Llon=Llon-dx;
          Rlon=max(rlon(:));  Rlon=Rlon+dx;
          Blat=min(rlat(:));  Blat=Blat-dy;
          Tlat=max(rlat(:));  Tlat=Tlat+dy;

          [C]=r_gshhs(Llon,Rlon,Blat,Tlat,Cname);
          [C]=x_gshhs(Llon,Rlon,Blat,Tlat,C,'patch');

          [y,x]=meshgrid(1:Mp,1:Lp);

          xcst=griddata(rlon,rlat,x,C.lon,C.lat,'linear');
          ycst=griddata(rlon,rlat,y,C.lon,C.lat,'linear');

          xcst=xcst-1;              % substract one to have indices
          ycst=ycst-1;              % in the range (0:L,0:M)

          save(ijcoast, 'xcst', 'ycst')
          disp(['Saved ij-coastline file:   ', ijcoast]);

          clear C x y
        else
          load(ijcoast);
          disp(['Loaded ij-coastline file:  ', ijcoast]);
        end
      elseif (got_coast)            % get coastline from GRID NetCDF file

        ijcoast=strcat(fname, '_jcst.mat');

        if (~exist(ijcoast, 'file'))
          C.lon=G.lon_coast;
          C.lat=G.lat_coast;

          [y,x]=meshgrid(1:Mp,1:Lp);

          xcst=griddata(rlon,rlat,x,C.lon,C.lat,'linear');
          ycst=griddata(rlon,rlat,y,C.lon,C.lat,'linear');

          xcst=xcst-1;              % substract one to have indices
          ycst=ycst-1;              % in the range (0:L,0:M)

          save(ijcoast, 'xcst', 'ycst')
          disp(['Saved ij-coastline file:   ', ijcoast]);

          clear C x y
        else
          load(ijcoast);
          disp(['loaded ij-coastline file:  ', ijcoast]);
        end

      elseif (got_clat && got_clon) % use provided coastline data

        ijcoast=strcat(fname, '_ijcst.mat');

        if (~exist(ijcoast, 'file'))
          [y,x]=meshgrid(1:Mp,1:Lp);

          xcst=griddata(rlon,rlat,x,C.lon,C.lat,'linear');
          ycst=griddata(rlon,rlat,y,C.lon,C.lat,'linear');

          xcst=xcst-1;              % substract one to have indices
          ycst=ycst-1;              % in the range (0:L,0:M)

          save(ijcoast, 'xcst', 'ycst')
          disp(['Saved ij-coastline file:   ', ijcoast]);

          clear x y
        else
          load(ijcoast);
          disp(['loaded ij-coastline file:  ', ijcoast]);
        end

      elseif (got_matfile)          % use provided .mat data

        load(coast_file);

        if (exist('C','var'))
          xcst=C.Icst;
          ycst=C.Jcst;
          clear C
        elseif (exist('lon','var') && exist('lat','var'))

          ijcoast=strcat(fname, '_ijcst.mat');          

          if (~exist(ijcoast, 'file'))
            [y,x]=meshgrid(1:Mp,1:Lp);

            xcst=griddata(rlon,rlat,x,lon,lat,'linear');
            ycst=griddata(rlon,rlat,y,lon,lat,'linear');

            xcst=xcst-1;            % substract one to have indices
            ycst=ycst-1;            % in the range (0:L,0:M)

            save(ijcoast, 'xcst', 'ycst')
            disp(['Saved ij-coastline file:   ', ijcoast]);

            clear lon lat x y
          else
            load(ijcoast);
            disp(['loaded ij-coastline file:  ', ijcoast]);
          end
        
        else
          error([coast_file, 'should contain "lon" and "lat" vectors']);
        end
      
      else
        error('unable to process coastline data');
      end

      toc;
      disp(' ');

    end

% Initialize the window.

    fig=figure('NumberTitle','off',                                     ...
               'tag','MaskEditor',                                      ...
               'DoubleBuffer','on',                                     ...
               'backingstore','off',                                    ...
               'menubar','none');

% Since starting with version 2013a Matlab no longer allows resizing
% figure by window manager, the following segment is to set optimal
% sizes by Matlab itself.  The working portion of the figure sized
% [wdth * hgth] which are set to be proportionally to the dimensions
% of the grid (thus always keeping the aspect ratio) and a specified
% portion of the screen size; after that left, right, bottom, and top
% side margins dL,dR,dB,dT are added to accommodate axes and menu
% buttons on the right side. These are specified in pixels and are
% kept constant regardless of figure size.

    dL=50; dR=220;   dB=45; dT=25;

    set(0,'Units','pixels')
    screen=get(0,'ScreenSize');

    hgth=0.9*screen(4);       wdth=hgth*Lp/Mp;
    full_width=wdth+dL+dR;    full_height=hgth+dB+dT;    

    xo=16;                                 % Working area left-bottom
    yo=screen(4)-full_height-40;           % corner coordinates, width,
    wb=full_width;                         % and height (in pixels)
    hb=full_height;
    set(fig,'Position', [xo yo wb hb]);

    [mx,my]=ndgrid(1:Lp,1:Mp);

    gx=mx(:,1)-1;
    gy=my(1,:)-1;    
    
    cdata=rmask*2+mod(mod(mx,2)+mod(my,2),2)+1;

    xn=dL/full_width;                      % Working area left-bottom
    yn=dB/full_height;                     % corner coordinates, width,
    wn=wdth/full_width;                    % and height (normalized)
    hn=hgth/full_height;
    axes('position', [xn yn wn hn]);
    hplot=image(gx,gy,cdata','cdatamapping','direct');
    set(gca,'YDir','normal',                                            ...
            'layer','top',                                              ...
            'tickdir','out');
    colormap(CMAP);
    hold on;
    hline=plot(xcst,ycst,LineColor);
    set(hline,'LineWidth',LineWidth);
    xl=xlim;
    yl=ylim;
    changed=0;
    setgui(full_width,wdth,dR,dL);
    editmask refresh;

end

return


function setgui(full_width,wdth,dR,dL)

%--------------------------------------------------------------------------
% set-up Land/Sea mask editing menu bottons.
%--------------------------------------------------------------------------

xpos=(dL+wdth+.05*dR)/full_width;
bwdth=0.9*dR/full_width;

textbox([xpos .85  bwdth .10], '(i,j)', {'0,0'}, 'pos');

radiobox([xpos .625 bwdth .19],                                         ...
         'Edit Mode',{'Land/Sea','Set Land','Set Sea'},'mode');

radiobox([xpos .45 bwdth .14],                                          ...
         'Edit Tool',{'Point edit','Area edit'},'tool');

button([xpos .375 bwdth .05],                                           ...
       'Zoom In',[mfilename ' zoomin']);

button([xpos .325 bwdth .05],                                           ...
       'Zoom Out',[mfilename ' zoomout']);

button([xpos .225 bwdth .05],                                           ...
       'Undo',[mfilename ' undo']);

button([xpos .175 bwdth .05],                                           ...
       'Save',[mfilename ' save']);

button([xpos .075 bwdth .05],                                           ...
       'Exit',[mfilename ' exit']);

axisscroll('r')
axisscroll('t')
set(gcf,'WindowButtonMotionFcn',[mfilename ' move;'],                   ...
        'CloseRequestFcn',[mfilename ' exit;'],                         ...
        'interruptible','on');
enable_click;

return

function disable_click

%--------------------------------------------------------------------------
% Disable pointer clicking on current figure.
%--------------------------------------------------------------------------

set(gcf,'WindowButtonDownFcn','');

return

function enable_click

%--------------------------------------------------------------------------
% Enable pointer clicking on current figure.
%--------------------------------------------------------------------------

set(gcf,'WindowButtonDownFcn',[mfilename ' click;']);

return

function r=within(a,b)

%--------------------------------------------------------------------------
% Check if 'a' is within the range of 'b'.
%--------------------------------------------------------------------------

r=(a>=b(1) & a<= b(2));

return
