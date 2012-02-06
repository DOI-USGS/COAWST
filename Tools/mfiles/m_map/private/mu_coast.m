function [ncst,Area,k]=mu_coast(optn,varargin);
% MU_COAST Add a coastline to a given map.
%         MU_COAST draw a coastline as either filled patches (slow) or
%         lines (fast) on a given projection. It uses a coastline database with
%         a resolution of about 1/4 degree. 
%
%         It is not meant to be called directly by the user, instead it contains
%         code common to various functions in the m_map directory.
%
%         Calling sequence is MU_COAST(option,arguments) where
%
%         Option string  
%           c,l,i,h,f :  Accesses various GSHHS databases. Next argument is
%                        GSHHS filename.
%
%           u[ser] :   Accesses user-specified coastline file (a mat-file of
%                      data saved from a previous MU_COAST call). Next argument
%                      is filename
%
%           v[ector] : Uses input vector of data. Next argument is the 
%                     data in the form of a nx2 matrix of [longitude latitude].
%                     Patches must have first and last points the same. In a vector,
%                     different patches can be separated by NaN.
%
%           d[efault] :  Accesses default coastline.
%
%         The arguments given above are then followed by (optional) arguments 
%         specifying lines or patches, in the form:
%
%          optional arguments:  <line property/value pairs>, or
%                               'line',<line property/value pairs>.
%                               'patch',<patch property/value pairs>.
%
%         If no or one output arguments are specified then the coastline is drawn, with
%         patch handles returned. 
%         If 3 output arguments are specified in the calling sequence, then no drawing
%         is done. This can be used to save subsampled coastlines for future use
%         with the 'u' option, for fast drawing of multiple instances of a particular
%         coastal region.
%
%
%    
%         See also M_PROJ, M_GRID     

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% Notes: 15/June/98 - fixed some obscure problems that occured sometimes
%                     when using conic projections with large extents that
%                     crossed the 180-deg longitude line.
%        31/Aug/98  - added "f" gshhs support (Thanks to RAMO)
%        17/June/99 - 'line' option as per manual (thanks to Brian Farrelly)
%        3/Aug/01   - kludge fix for lines to South Pole in Antarctica (response
%                     to Matt King).
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.


global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;


% m_coasts.mat contains 3 variables:
% ncst: a Nx2 matrix of [LONG LAT] line segments, each of which form
%       a closed contour, separated by NaN
%  k=[find(isnan(ncst(:,1)))];
% Area: a vector giving the areas of the different patches. Both ncst
%     and Area should be ordered with biggest regions first. Area can
%     be computed as follows:
%
%  Area=zeros(length(k)-1,1);
%  for i=1:length(k)-1,
%    x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
%    y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
%    nl=length(x);
%    Area(i)=sum( diff(x).*(y(1:nl-1)+y(2:nl))/2 );
%  end;
%
%     Area should be >0 for land, and <0 for lakes and inland seas.

switch optn(1),
  case {'c','l','i','h','f'},  
    [ncst,k,Area]=get_coasts(optn,varargin{1});
    varargin(1)=[];
  case  'u',                   
    eval(['load ' varargin{1} ' -mat']);
    varargin(1)=[];
  case 'v',
    ncst=[NaN NaN;varargin{1};NaN NaN];
    varargin(1)=[];
    k=[find(isnan(ncst(:,1)))];  % Get k
    Area=ones(size(k));          % Make dummy Area vector (all positive).
  otherwise
    load m_coasts
end;

% If all we wanted was to extract a sub-coastline, return.
if nargout==3,
  return;
end;

% Handle wrap-arounds (not needed for azimuthal and oblique projections)

switch MAP_PROJECTION.routine,
 case {'mp_cyl','mp_conic','mp_tmerc'}
  if MAP_VAR_LIST.longs(2)<-180,
   ncst(:,1)=ncst(:,1)-360;
  elseif MAP_VAR_LIST.longs(1)>180,
   ncst(:,1)=ncst(:,1)+360;
  elseif MAP_VAR_LIST.longs(1)<-180,
   Area=[Area;Area];
   k=[k;k(2:end)+k(end)-1];
   ncst=[ncst;[ncst(2:end,1)-360 ncst(2:end,2)]];
   % This is kinda kludgey - but sometimes adding all these extra points
   % causes wrap-around in the conic projection, so we want to limit the
   % longitudes to the range needed. However, we don't just clip them to
   % min long because that can cause problems in trying to decide which way
   % curves are oriented when doing the fill algorithm below. So instead
   % I sort of crunch the scale, preserving topology.
   nn=ncst(:,1)<MAP_VAR_LIST.longs(1);
   ncst(nn,1)=(ncst(nn,1)-MAP_VAR_LIST.longs(1))/100+MAP_VAR_LIST.longs(1);
  elseif MAP_VAR_LIST.longs(2)>180,
   Area=[Area;Area];
   k=[k;k(2:end)+k(end)-1];
   ncst=[ncst;[ncst(2:end,1)+360 ncst(2:end,2)]];
   % Ditto.
   nn=ncst(:,1)>MAP_VAR_LIST.longs(2);
   ncst(nn,1)=(ncst(nn,1)-MAP_VAR_LIST.longs(2))/100+MAP_VAR_LIST.longs(2);
  end;
end;


if length(varargin)>0,
  if strcmp(varargin(1),'patch'),
    optn='patch';
  end
  if strcmp(varargin(1),'line'),
    optn='line';
    varargin=varargin(2:end); % ensure 'line' does not get passed to line
  end
else
  optn='line';
end;



switch optn,
 case 'patch',

  switch MAP_VAR_LIST.rectbox,
    case 'on',
      xl=MAP_VAR_LIST.xlims;
      yl=MAP_VAR_LIST.ylims;
      [X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','on');
      oncearound=4;
    case 'off',
      xl=MAP_VAR_LIST.longs;
      yl=MAP_VAR_LIST.lats;
      X=ncst(:,1);
      Y=ncst(:,2);
      [X,Y]=mu_util('clip','on',X,xl(1),X<xl(1),Y);
      [X,Y]=mu_util('clip','on',X,xl(2),X>xl(2),Y);
      [Y,X]=mu_util('clip','on',Y,yl(1),Y<yl(1),X);
      [Y,X]=mu_util('clip','on',Y,yl(2),Y>yl(2),X);
      oncearound=4;
    case 'circle',
      rl=MAP_VAR_LIST.rhomax;
      [X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','on');
      oncearound=2*pi;    
  end;

  p_hand=zeros(length(k)-1,1); % Patch handles
  
  for i=1:length(k)-1,
    x=X(k(i)+1:k(i+1)-1);
    fk=finite(x);
    if any(fk),
      y=Y(k(i)+1:k(i+1)-1);
      nx=length(x);
      if Area(i)<0, x=flipud(x);y=flipud(y); fk=flipud(fk); end;
%clf
%line(x,y,'color','m');
      st=find(diff(fk)==1)+1;
      ed=find(diff(fk)==-1);
%length(x),
%ed,
%st
       if length(st)<length(ed) | isempty(st), st=[ 1;st]; end;
       if length(ed)<length(st),               ed=[ed;nx]; end;
%ed
%st
      if  ed(1)<st(1),
        if c_edge(x(1),y(1))==9999,
          x=x([(ed(1)+1:end) (1:ed(1))]);
          y=y([(ed(1)+1:end) (1:ed(1))]);
          fk=finite(x);
          st=find(diff(fk)==1)+1;
          ed=[find(diff(fk)==-1);nx];
          if length(st)<length(ed), st=[1;st]; end
        else
          ed=[ed;nx];
          st=[1;st];
        end;
      end;
%ed
%st
      % Get rid of 2-point curves (since often these are out-of-range lines with
      % both endpoints outside the region of interest)
      k2=(ed-st)<3;
      ed(k2)=[]; st(k2)=[];
%%[XX,YY]=m_ll2xy(x(st),y(st),'clip','off');
%line(x,y,'color','r','linest','-');
%line(x(st),y(st),'marker','o','color','r','linest','none');
%line(x(ed),y(ed),'marker','o','color','g','linest','none');
      edge1=c_edge(x(st),y(st));
      edge2=c_edge(x(ed),y(ed));
%-edge1*180/pi
%-edge2*180/pi
      mi=1;
      while length(st)>0,
        if mi==1,
          xx=x(st(1):ed(1));
          yy=y(st(1):ed(1));
        end;
        estart=edge2(1);
        s_edge=edge1;
        s_edge(s_edge<estart)=s_edge(s_edge<estart)+oncearound;
%s_edge,estart
        [md,mi]=min(s_edge-estart);
        switch MAP_VAR_LIST.rectbox,
          case {'on','off'},
            for e=floor(estart):floor(s_edge(mi)),
              if e==floor(s_edge(mi)), xe=x(st([mi mi])); ye=y(st([mi mi])); 
              else  xe=xl; ye=yl; end;
              switch rem(e,4),
                case 0,
                  xx=[xx; xx(end*ones(10,1))];
                  yy=[yy; yy(end)+(ye(2)-yy(end))*[1:10]'/10 ];
                case 1,
                  xx=[xx; xx(end)+(xe(2)-xx(end))*[1:10]'/10 ];
                  yy=[yy; yy(end*ones(10,1))];
                case 2,
                  xx=[xx; xx(end*ones(10,1))];
                  yy=[yy; yy(end)+(ye(1)-yy(end))*[1:10]'/10;];
                case 3,
                  xx=[xx; xx(end)+(xe(1)-xx(end))*[1:10]'/10 ];
                  yy=[yy; yy(end*ones(10,1))];
              end;
            end;
          case 'circle',
            if estart~=9999,
%s_edge(mi),estart
              xx=[xx; rl*cos(-[(estart:.1:s_edge(mi))]')];
              yy=[yy; rl*sin(-[(estart:.1:s_edge(mi))]')];
            end;
        end;
        if mi==1, % joined back to start
          if strcmp(MAP_VAR_LIST.rectbox,'off'), 
            [xx,yy]=m_ll2xy(xx,yy,'clip','off'); 
          end;
%disp(['paused-1 ' int2str(i)]);pause;
          if Area(i)<0,
            p_hand(i)=patch(xx,yy,varargin{2:end},'facecolor',get(gca,'color'));
          else
            p_hand(i)=patch(xx,yy,varargin{2:end});
          end;
          ed(1)=[];st(1)=[];edge2(1)=[];edge1(1)=[];
        else
          xx=[xx;x(st(mi):ed(mi))];
          yy=[yy;y(st(mi):ed(mi))];
          ed(1)=ed(mi);st(mi)=[];ed(mi)=[];
          edge2(1)=edge2(mi);edge2(mi)=[];edge1(mi)=[];
        end;
%%disp(['paused-2 ' int2str(i)]);pause;
      end;
    end; 
  end;  

 otherwise,
 
  % This handles the odd points required at the south pole by any Antarctic
  % coastline by setting them to NaN (for lines only)
  ii=ncst(:,2)<=-89.9;
  if any(ii), ncst(ii,:)=NaN; end;
  
  [X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','on');
 
  % Get rid of 2-point lines (these are probably clipped lines spanning the window)
  fk=finite(X);        
  st=find(diff(fk)==1)+1;
  ed=find(diff(fk)==-1);
  k=find((ed-st)==1);
  X(st(k))=NaN;

  p_hand=line(X,Y,varargin{:}); 

end;

ncst=p_hand;

%-----------------------------------------------------------------------
function edg=c_edge(x,y);
% C_EDGE tests if a point is on the edge or not. If it is, it is
%        assigned a value representing it's position oon the perimeter
%        in the clockwise direction. For x/y or lat/long boxes, these
%        values are
%           0 -> 1 on left edge
%           1 -> 2 on top
%           2 -> 3 on right edge
%           3 -> 4 on bottom
%        For circular boxes, these values are the -ve of the angle
%        from center.

global MAP_VAR_LIST

switch MAP_VAR_LIST.rectbox,
  case 'on',
    xl=MAP_VAR_LIST.xlims;
    yl=MAP_VAR_LIST.ylims;
  case 'off',
    xl=MAP_VAR_LIST.longs;
    yl=MAP_VAR_LIST.lats;
  case 'circle',
    rl2=MAP_VAR_LIST.rhomax^2;
end;

edg=9999+zeros(length(x),1);
tol=1e-10;

switch MAP_VAR_LIST.rectbox,
  case {'on','off'},
    i=abs(x-xl(1))<tol;
    edg(i)=(y(i)-yl(1))/diff(yl);

    i=abs(x-xl(2))<tol;
    edg(i)=2+(yl(2)-y(i))/diff(yl);

    i=abs(y-yl(1))<tol;
    edg(i)=3+(xl(2)-x(i))/diff(xl);

    i=abs(y-yl(2))<tol;
    edg(i)=1+(x(i)-xl(1))/diff(xl);

  case 'circle',
    i=abs(x.^2 + y.^2 - rl2)<tol;
    edg(i)=-atan2(y(i),x(i));   % use -1*angle so that numeric values
                                % increase in CW direction
end;


%%
function [ncst,k,Area]=get_coasts(optn,file);
%
%  GET_COASTS  Loads various GSHHS coastline databases and does some preliminary
%              processing to get things into the form desired by the patch-filling
%              algorithm.
%
% Changes; 3/Sep/98 - RP: better decimation, fixed bug in limit checking.

global MAP_PROJECTION MAP_VAR_LIST

llim=rem(MAP_VAR_LIST.longs(1)+360,360)*1e6;
rlim=rem(MAP_VAR_LIST.longs(2)+360,360)*1e6;
tlim=MAP_VAR_LIST.lats(2)*1e6;
blim=MAP_VAR_LIST.lats(1)*1e6;

mrlim=rem(MAP_VAR_LIST.longs(2)+360+180,360)-180;
mllim=rem(MAP_VAR_LIST.longs(1)+360+180,360)-180;
mtlim=MAP_VAR_LIST.lats(2);
mblim=MAP_VAR_LIST.lats(1);

% decfac is for decimation of areas outside the lat/long bdys.
switch optn(1),
  case 'f',   % 'full' (undecimated) database
    ncst=NaN+zeros(492283,2);Area=zeros(41520,1);k=ones(41521,1);
    decfac=12500;
  case 'h',
    ncst=NaN+zeros(492283,2);Area=zeros(41520,1);k=ones(41521,1);
    decfac=2500;
  case 'i',
    ncst=NaN+zeros(492283,2);Area=zeros(41520,1);k=ones(41521,1);
    decfac=500;
  case 'l',
    ncst=NaN+zeros(101023,2);Area=zeros(10768,1);k=ones(10769,1);
    decfac=100;
  case 'c',
    ncst=NaN+zeros(14872,2);Area=zeros(1868,1);k=ones(1869,1);
    decfac=20;
end;
fid=fopen(file,'r','ieee-be');

if fid==-1,
  warning(sprintf(['Coastline file ' file ...
          ' not found \n(Have you installed it? See the M_Map User''s Guide for details)' ...
          '\n ---Using default coastline instead']));
  load m_coasts
  return
end;


Area2=Area;

[A,cnt]=fread(fid,9,'int32');

l=0;
while cnt>0,

 % A: 1:ID, 2:num points, 3:land/lake etc., 4-7:w/e/s/n, 8:area, 9:greenwich crossed.
 
 C=fread(fid,A(2)*2,'int32'); % Read all points in the current segment.

 
 a=rlim>llim;  % Map limits cross longitude jump? (a==1 is no)
 b=A(9)<65536; % Cross boundary? (b==1 if no).
 c=llim<rem(A(5)+360e6,360e6); 
 d=rlim>rem(A(4)+360e6,360e6);
 e=tlim>A(6) & blim<A(7);
 
 % This test checks whether the lat/long box containing the line overlaps that of
 % the map. There are various cases to consider, depending on whether map limits
 % and/or the line limits cross the longitude jump or not.
 
%% if e & (  a&( b&c&d | ~b&(c|d)) | ~a&(~b | (b&(c|d))) ),
 if e & (  (a&( (b&c&d) | (~b&(c|d)) )) | (~a&(~b | (b&(c|d))) ) ),
 
   l=l+1;
 
   x=C(1:2:end)*1e-6;y=C(2:2:end)*1e-6;

   %  make things continuous (join edges that cut across 0-meridian)

   dx=diff(x);
%%fprintf('%f %f\n', max(dx),min(dx))
%   if A(9)>65536 | any(dx)>356 | any(dx<356),
%  if ~b | any(dx>356) | any(dx<-356),
     x=x-360*cumsum([x(1)>180;(dx>356) - (dx<-356)]);
%   end;

   % Antarctic is a special case - extend contour to make nice closed polygon
   % that doesn't surround the pole.   
   if abs(x(1))<1 & abs(y(1)+68.9)<1,
     y=[-89.9;-78.4;y(x<=-180);y(x>-180);   -78.4;-89.9*ones(18,1)];
     x=[  180; 180 ;x(x<=-180)+360;x(x>-180);-180; [-180:20:160]'];
   end;

   % First and last point should be the same.
   
   if x(end)~=x(1) | y(end)~=y(1), x=[x;x(1)];y=[y;y(1)]; end;

   % get correct curve orientation for patch-fill algorithm.
   
   Area2(l)=sum( diff(x).*(y(1:(end-1))+y(2:end))/2 );
   Area(l)=A(8)/10;

   if rem(A(3),2)==0; 
     Area(l)=-abs(Area(l)); 
     if Area2(l)>0, x=x(end:-1:1);y=y(end:-1:1); end;
   else
     if Area2(l)<0, x=x(end:-1:1);y=y(end:-1:1); end; 
   end;

   % Here we try to reduce the number of points.
   
   xflag=0;
   if max(x)>180, % First, save original curve for later if we anticipate
     sx=x;sy=y;   % a 180-problem.
     xflag=1;
   end;
   
   % Look for points outside the lat/long boundaries, and then decimate them
   % by a factor of about 'decfac' (don't get rid of them completely because that
   % can sometimes cause problems when polygon edges cross curved map edges).
   
   tol=.2;   
  
   % Do y limits, then x so we can keep corner points.
   
   nn=y>mtlim+tol | y<mblim-tol;
     % keep one extra point when crossing limits, also the beginning/end point.
   nn=logical(nn-([0;diff(nn)]>0)-([diff(nn);0]<0));nn([1 end])=0;
     % decimate vigorously
   nn=nn & rem(1:length(nn),decfac)'~=0;
   x(nn)=[];y(nn)=[];
         
   if mrlim>mllim,  % no wraparound
       % sections of line outside lat/long limits
     nn=(x>mrlim+tol | x<mllim-tol) & y<mtlim & y>mblim;
    else            % wraparound case
     nn=(x>mrlim+tol & x<mllim-tol ) & y<mtlim & y>mblim;
   end;
   nn=logical(nn-([0;diff(nn)]>0)-([diff(nn);0]<0));nn([1 end])=0;
   nn=nn & rem(1:length(nn),decfac)'~=0;
   x(nn)=[];y(nn)=[];
   
   % Move all points "near" to map boundaries.
   % I'm not sure about the wisdom of this - it might be better to clip
   % to the boundaries instead of moving. Hmmm. 
      
   y(y>mtlim+tol)=mtlim+tol;
   y(y<mblim-tol)=mblim-tol;
   if mrlim>mllim,   % Only clip long bdys if I can tell I'm on the right
                     % or left (i.e. not in wraparound case)
     x(x>mrlim+tol)=mrlim+tol;
     x(x<mllim-tol)=mllim-tol;
   end;   
   
 %% plot(x,y);pause;
 
   k(l+1)=k(l)+length(x)+1;
   ncst(k(l)+1:k(l+1)-1,:)=[x,y];
   ncst(k(l+1),:)=[NaN NaN];
   
   % This is a little tricky...the filling algorithm expects data to be in the
   % range -180 to 180 deg long. However, there are some land parts that cut across
   % this divide so they appear at +190 but not -170. This causes problems later...
   % so as a kludge I replicate some of the problematic features at 190-360=-170.
   % Small islands are just duplicated, for the Eurasian landmass I just clip
   % off the eastern part.
   
   if xflag,
     l=l+1;Area(l)=Area(l-1); 
     if abs(Area(l))>1e5,
       nn=find(sx>180);nn=[nn;nn(1)];
       k(l+1)=k(l)+length(nn)+1;
       ncst(k(l)+1:k(l+1)-1,:)=[sx(nn)-360,sy(nn)];
     else   % repeat the island at the other edge.
       k(l+1)=k(l)+length(sx)+1;
       ncst(k(l)+1:k(l+1)-1,:)=[sx-360,sy];
     end;
     ncst(k(l+1),:)=[NaN NaN];
   end;
 end;
 
 
 [A,cnt]=fread(fid,9,'int32');

end;

fclose(fid);

%%plot(ncst(:,1),ncst(:,2));pause;clf;
%size(ncst)
%size(Area)
%size(k)


ncst((k(l+1)+1):end,:)=[];  % get rid of unused part of data matrices
Area((l+1):end)=[];
k((l+2):end)=[];































