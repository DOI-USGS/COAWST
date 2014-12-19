function editbath(grid_file, coast_file);

% EDITMASK:  Interactive Land/Sea mask editing for ROMS
%
% EDITMASK(GRID_FILE,COAST_FILE) is tool for manual editing of the
% Land/Sea mask on RHO-points.  To accelerate the proccessing, the
% Land/Sea mask is edited in (I,J) grid coordinates.  GRID_FILE is
% the GRID NetCDF file containing the grid and mask. COAST_FILE is
% a MAT file holding the coastline (lon,lat) coordinates  or (I,J)
% grid coordinates.  If the (I,J) coordinates are not provided, it
% will compute and write them into file. If called without, one or
% both arguments, it will prompt for the needed file name(s).
%
% Mouse shortcuts:
%
%    double click ==> Zoom in
%    right  click ==> Zoom out
%    middle click ==> change editing mode
%
%    Calls: READ_MASK, WRITE_MASK and UVP_MASKS functions.
%           BUTTON, RADIOBOX, TEXTBOX, AXISSCROLL,
%

% svn $Id: editbath.m 436 2010-01-02 17:07:55Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                            A. Shcherbina          %
%===========================================================================%

% Define and initialize persistent variables.

persistent changed rmask rbath rlon rlat bath mask hplot hcst ha Lp Mp mx my
persistent mfile fig zooming xx yy xl yl xcst ycst

% Single global variable to pass info to/from the callback routines.

global GUI

% Set colormap: first two entries - land, second pair - sea
CMAP=[.5 1 0;1 1 0;0 0 .7;0 0 1];

% Set coastline line color and width.

LineColor='k';
LineWidth=1;

% Check input arguments.

if (nargin < 1 || isempty(grid_file))
  grid_file='*.nc';
end
if (nargin < 2 || isempty(coast_file))
  coast_file='*.mat';
  got_coast=1;
end

FIGURE_NAME='Land/Sea and Bathymetry Editor';

%===========================================================================
% The EDITMASK function is also used as a CallBack for some of the
% uicontrols,  so we need to figure out, what's required by looking
% at the 'grid_file' parameter.
%===========================================================================

switch lower(grid_file)

%---------------------------------------------------------------------------
% Zoom-in.
%---------------------------------------------------------------------------

  case 'zoomin'

    disable_click;
    editbath move;
    zooming=1;
    waitforbuttonpress
    xx0=xx; yy0=yy;                             % save current pointer position
    rbbox;                                      % track rubberband rectangle
    xx1=xx; yy1=yy;                             % new pointer position
    if ((xx0 ~= xx1) && (yy0 ~= yy1))           % trim limits and set zoom
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

%---------------------------------------------------------------------------
% Zoom-out.
%---------------------------------------------------------------------------

  case 'zoomout'

    xlim(xl);
    ylim(yl);
    axisscroll;

%---------------------------------------------------------------------------
% Edit Land/Sea and bathymetry.
%---------------------------------------------------------------------------

  case 'click'

    button=get(gcf, 'SelectionType');
    if (strcmp(button,'alt'))                  % zoom out on right click
      editbath zoomout;
      return;
    end
    if (strcmp(button,'open'))                 % zoom in on double click
      editbath zoomin;
      return;
    end,
    if (strcmp(button,'extend'))               % cycle modes on middle click
      m=mod(GUI.mode,3)+1;
      eval(get(GUI.mode_h(m),'callback'));
      editbath zoomout;
      return;
    end
    if (within(xx,xlim) && within(yy,ylim))     % left click within edit area
      disable_click;
      switch GUI.tool
        case 1                                  % point edit
          ix=floor(xx+1.5);
          iy=floor(yy+1.5);
          switch GUI.mode
            case 1                             % toggle between land and sea
              rmask(ix,iy)=~rmask(ix,iy);
            case 2                             % set bathymetry
              if rmask(ix,iy)==1
              prompt = {'Enter bathymetry: '};
              num_lines = 1;
              dlg_title = ['Modify point: ',...
                            num2str(ix),' ',num2str(iy)];
              answer = inputdlg(prompt,dlg_title,num_lines);
              rbath(ix,iy)=str2double(answer);
              else
              rbath(ix,iy)= 0;    
              end
          end
        case 2                                  % area edit
          xx0=xx; yy0=yy;                       % save current pointer position
          rbbox;                                % track rubberband rectangle
          xx1=xx; yy1=yy;                       % new pointer position
          idsel=find((mx-xx0-1).*(mx-xx1-1)<=0 & ...
                     (my-yy0-1).*(my-yy1-1)<=0);% indicies within rectangle
          switch GUI.mode
            case 1                              % toggle between land and sea
              rmask(idsel)=~rmask(idsel);       
            case 2                              % set bathymetry
              prompt = {'Enter bathymetry: '};
              num_lines = 1;
              dlg_title = 'Modify area ';
              answer = inputdlg(prompt,dlg_title,num_lines);
              rbath(idsel)=str2double(answer);
          end
      end
      changed=1;
      enable_click;
      editbath refresh;
    end

%---------------------------------------------------------------------------
% Update mask by changing colormap.
%---------------------------------------------------------------------------

  case 'refresh',

    cdata=rmask*2+mod(mod(mx,2)+mod(my,2),2)+1;
    set(hplot,'cdata',cdata');
    nm=[FIGURE_NAME,' - ', mfile];
    if (changed),
      nm=[nm,'*'];
    end,
    set(fig,'name',nm);

%---------------------------------------------------------------------------
% Pointer movement: update pointer coordinates on menu.
%---------------------------------------------------------------------------

  case 'move',

    xy=get(gca,'currentpoint');
    xx=xy(1,1); yy=xy(1,2);
    if (within(xx,xlim) && within(yy,ylim))
      bb=rbath(floor(xx+1.5),floor(yy+1.5));
      s=sprintf('[%d,%d,%f]',floor(xx+.5),floor(yy+.5),bb);
      if (zooming),
        pointer('zoom');
      elseif (GUI.tool==1),
        pointer('point');
      elseif GUI.tool==2,
        pointer('rect');
      end,
    else
      s='---';
      pointer('arrow');
    end
    set(GUI.pos_h,'string',s);

%---------------------------------------------------------------------------
% Compute U-, V- and PSI masks.  Write out into GRID NetCDF file.
%---------------------------------------------------------------------------

  case 'save',
    [umask,vmask,pmask]=uvp_masks(rmask);
    mask=rmask;
    bath=rbath;
    write_mask(mfile,rmask,umask,vmask,pmask);
    write_bath(mfile,rbath);
    changed=0;
    editbath refresh;

%---------------------------------------------------------------------------
% Undo changes: restore last saved mask.
%---------------------------------------------------------------------------

  case 'undo',
    rmask=mask;
    rbath=bath;
    changed=0;
    editbath refresh;

%---------------------------------------------------------------------------
% Done: inquire to save mask changes.
%---------------------------------------------------------------------------

  case 'exit',

    if (~changed),
      delete(gcf);
    else
      res=questdlg('The bathy has been changed. Save?',FIGURE_NAME);
      switch res
        case 'Yes',
          editbath save;
          disp('Mask and bathy have been saved');
          delete(gcf);
        case 'No',
          disp('Mask has NOT been saved');
          delete(gcf);
      end,
    end,

%---------------------------------------------------------------------------
% Help.
%---------------------------------------------------------------------------

  case 'help',

    show_help;

%---------------------------------------------------------------------------
% Initialize: read mask and coastline data.
%---------------------------------------------------------------------------

  otherwise,

% Kill all windows.

    delete(findobj(0,'tag','MaskEditor'));

% If appropriate, inquire input file names.

    if (any(grid_file=='*')),
      [fn,pth]=uigetfile(grid_file,'Select ROMS grid file...');
      if (~fn),
        return,
      end;
      grid_file=[pth,fn];
    end,
    if (any(coast_file=='*')),
      [fn,pth]=uigetfile(coast_file,'Select the coastline ...');
      if (~fn),
        return,
      end;
      coast_file=[pth,fn];
    end,
    

% Read in grid data.

   mfile=grid_file;
   [spherical,rlon,rlat,bath,mask]=read_mask(grid_file);
   rmask=mask;
   rbath=bath;
   [Lp,Mp]=size(mask);
   [mx,my]=ndgrid(1:Lp,1:Mp);

% Check if grid NetCDF file has coastline data.

   [varnam,nvars]=nc_vname(grid_file);
   got_Clon=strmatch('lon_coast',varnam,'exact');
   got_Clat=strmatch('lat_coast',varnam,'exact');
   if isempty(got_Clon),
     got_Clon=0;
   end
   if isempty(got_Clat),
     got_Clat=0;
   end

   if (~(got_Clon & got_Clat) & spherical),
     if (any(coast_file=='*')),
       [fn,pth]=uigetfile(coast_file,'Select Coastline file...');
       if (~fn),
         return,
         got_coast=0;
       end;
       coast_file=[pth,fn];
     end,
   end,

% Read in coastline data. If appropriate, compute coastline (I,J)
% grid indices.

   if (got_Clon && got_Clat)
     [C]=ijcoast(grid_file);
     xcst=C.Icst;
     ycst=C.Jcst;
     clear C lat lon;
   else
     if (got_coast)
       load(coast_file);
       if (exist('C'))
         xcst=C.Icst;
         ycst=C.Jcst;
         clear C;
       elseif (exist('lon') && exist('lat'))
         prompt = {'yes=1/no=0'};  
         defaultanswer={'1'};
         num_lines = 1;
         dlg_title = 'Do you want to regenerate the mask with the coastline ?';
         options.Resize='on';
         options.WindowStyle='normal';
         options.Interpreter='tex';
         answer = inputdlg(prompt,dlg_title,num_lines,defaultanswer,options);
         yesno= str2double(answer(1,:));
         if(yesno==1)
            rmask_temp=0;
            lon_rho=nc_read(grid_file,'lon_rho');
            lat_rho=nc_read(grid_file,'lat_rho');
            test_nan=find(~isnan(lat));
            test_close=find(test_nan(2:end)-test_nan(1:end-1)~=1);
            [ntestnan,~]=size(test_nan);
            test_close=[-1;test_close;ntestnan];
            [ntestclose,~]=size(test_close);
            for i=1:ntestclose-1
                Xtemp=lon(test_nan(test_close(i)+2):test_nan(test_close(i+1)));
                Ytemp=lat(test_nan(test_close(i)+2):test_nan(test_close(i+1))); 
                if(isempty(find(Xtemp(2:end)-Xtemp(1:end-1)~=0))||isempty(find(Ytemp(2:end)-Ytemp(1:end-1)~=0))) 
                    disp('This is not a polygon but a line ...')
                else
                    rmask_temp=rmask_temp+inpoly_V1([lon_rho(:) lat_rho(:)],[Xtemp Ytemp]);
                end
            end
            rmask_temp=-(rmask_temp-1);
            rmask=reshape(rmask_temp,Lp,Mp);
         end     
         [C]=ijcoast(grid_file,coast_file);
         xcst=C.Icst;
         ycst=C.Jcst;
         clear C lat lon;
       else
         error('Coast file should contain "lon" and "lat" vectors');
       end
     end
   end
   
   % Initialize the window.
   fig=figure('NumberTitle','off',...
              'tag','MaskEditor',...
              'DoubleBuffer','on',...
              'backingstore','off',...
              'menubar','none');
   cdata=rmask*2+mod(mod(mx,2)+mod(my,2),2)+1;
   gx=mx(:,1)-1;
   gy=my(1,:)-1;
   axes('position',[0.0554 0.1024 0.7518 0.8405]);
   hplot=image(gx,gy,cdata','cdatamapping','direct','erasemode','normal');
   set(gca,'YDir','normal',...
           'drawmode','fast',...
           'layer','top',...
           'tickdir','out');
   colormap(CMAP);
   hold on;
   hline=plot(xcst,ycst,LineColor);
   set(hline,'LineWidth',LineWidth);
   xl=xlim; yl=ylim;
   changed=0;
   setgui;
   editbath refresh;

end
return


function setgui
%---------------------------------------------------------------------------
% set-up Land/Sea mask editing menu bottons.
%---------------------------------------------------------------------------

textbox([.85 .85 .145 .14], ....
        '(i,j,z)',{'0,0,0'},'pos');

radiobox([.85 .65 .145 .19], ...
         'Edit Mode',{'Toggle Land/Sea','Set Bathy'},'mode');

radiobox([.85 .5 .145 .14], ...
         'Edit Tool',{'Point edit','Area edit'},'tool');

button([.85 .4 .145 .05], ...
       'Zoom In',[mfilename ' zoomin']);

button([.85 .35 .145 .05], ...
       'Zoom Out',[mfilename ' zoomout']);

button([.85 .25 .145 .05], ...
       'Undo',[mfilename ' undo']);

button([.85 .2 .145 .05], ...
       'Save',[mfilename ' save']);

button([.85 .1 .145 .05], ...
       'Exit',[mfilename ' exit']);

axisscroll('r')
axisscroll('t')
set(gcf,'WindowButtonMotionFcn',[mfilename ' move;'],...
        'CloseRequestFcn',[mfilename ' exit;'],...
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
%---------------------------------------------------------------------------
% Enable pointer clicking on current figure.
%---------------------------------------------------------------------------

set(gcf,'WindowButtonDownFcn',[mfilename ' click;']);
return

function r=within(a,b)
%---------------------------------------------------------------------------
% Check if 'a' is within the range of 'b'.
%---------------------------------------------------------------------------

r=(a>=b(1) & a<= b(2));

return
