function [C]=r_gshhs(Llon, Rlon, Blat, Tlat, Fname)

%
% R_GSHHS:  Read and extract coastline data from GSHHS dataset
%
% [C]=r_gshhs(Llon, Rlon, Blat, Tlat, Fname)
%
%  This function reads requested GSHHS coastline database and extracts
%  data within requested map corners.  It also does some preliminary
%  processing to get things into the form desired by the patch-filling
%  algorithm.
%
%  Adapted from Rich Pawlowicz "get_coasts" function.
%
%  On Input:
%
%     Llon      Left   corner longitude (West values are negative)
%
%     Rlon      Right  corner longitude (West values are negative)
%
%     Blat      Bottom corner latitude (South values are negative)
%
%     Tlat      Top    corner latitude (South values are negative)
%
%     Fname     GSHHS database file name (string).
%
%  On Output:
%
%     C         Read coastline data (structure array):
%                 C.lon  => longitude of closed polygons separated by NaNs
%                 C.lat  => latitude  of closed polygons separated by NaNs
%                 C.area => polygon areas
%                 C.type => polygon type:
%                           C.type=1 => land,
%                           C.type=2 => lake,
%                           C.type=3 => island_in_lake,
%                           C.type=4 => pond_in_island_in_lake
%                 C.k    => number of points in polygon
%

% svn $Id: r_gshhs.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  Set extracting coordinates.
%--------------------------------------------------------------------------

llim=rem(Llon+360,360)*1e6;          %   micro-degrees
rlim=rem(Rlon+360,360)*1e6;
blim=Blat*1e6;
tlim=Tlat*1e6;

m.llim=rem(Llon+360+180,360)-180;
m.rlim=rem(Rlon+360+180,360)-180;
m.blim=Blat;
m.tlim=Tlat;

%--------------------------------------------------------------------------
%  Initialize output data: "decfac" is for decimation of areas outside the
%  lat/lon boundarys.
%--------------------------------------------------------------------------

if (strfind(Fname,'gshhs_f'))          % Full resolution database
  ncst=NaN(10810000,2);
  Area=zeros(188611,1);
  type=ones(188612,1);
  k=ones(188612,1);
  decfac=12500;
elseif (strfind(Fname,'gshhs_h'))      % High resolution database
  ncst=NaN(2080000,2);
  Area=zeros(153545,1);
  type=ones(153546,1);
  k=ones(153546,1);
  decfac=2500;
elseif (strfind(Fname,'gshhs_i'))      % Intermediate resolution database
  ncst=NaN(493096,2);
  Area=zeros(41529,1);
  type=ones(41530,1);
  k=ones(41530,1);
  decfac=500;
elseif (strfind(Fname,'gshhs_l'))      % Low resolution database
  ncst=NaN(124871,2);
  Area=zeros(20776,1);
  type=ones(27524,1);
  k=ones(27524,1);
  decfac=100;
elseif (strfind(Fname,'gshhs_c'))      % Crude resolution databse
  ncst=NaN(14403,2);
  Area=zeros(1766,1);
  type=ones(1767,1);
  k=ones(1767,1);
  decfac=20;
else
  error(['READ_GSHHS - unable to process database: ',Fname]);
end

flaglim=9;

%--------------------------------------------------------------------------
% Read in GSHHS database.
%--------------------------------------------------------------------------

fid=fopen(Fname,'r','ieee-be');

if (fid == -1)
  error(['R_GSHHS - unable to open file: ', Fname]);
end

Area2=Area;

% Read in header (A) for first poligon:
%
%      A(1):   polygn ID number
%      A(2):   Number of points
%      A(3):   1=land, 2=lake, 3=island_in_lake, 4=pond_in_island_in_lake
%      A(4):   Minimum longitude, w (micro-degrees)
%      A(5):   Maximum longitude, e (micro-degrees)
%      A(6):   Minimum latitude, s (micro-degrees)
%      A(7):   Maximum latitude, n (micro-degrees)
%      A(8):   Area of polygon (1/10 km^2)
%      A(9):   1 if greenwich is crossed

[cnt,g]=get_gheader(fid);

l=0;

while (cnt > 0)

% Read all points in the current segment. 

  C=fread(fid,g.N*2,'int32');    % Read all points in the current segment.
  
% For versions > 12 there are 2 Antarctics - ice-line and grounding line,
% also they changed the limits from 0-360 to -180 to 180.

  if (g.ver > 14)
    switch g.level
      case 6   
        g.level=flaglim+1; % ice line - don't show
      case 5   
        g.level=1;         % grounding line - use this.
    end

% For Antarctica the lime limits are 0 to 360 (exactly), thus c==0 and the
% line is not chosen for (e.g. a conic projection of part of Antarctica)
  
  else
    if (g.extentE == 360e6)
      g.extentE=g.extentE-1;
    end 
  end
  
  a=rlim > llim;            % Map limits cross longitude jump? (a==1 is no)
  b=g.greenwich < 65536;    % Cross boundary? (b==1 if no).
  c=llim < rem(g.extentE+360e6,360e6); 
  d=rlim > rem(g.extentW+360e6,360e6);
  e=tlim > g.extentS & blim < g.extentN;  

% This test checks whether the lat/long box containing the line overlaps
% that of the map. There are various cases to consider, depending on
% whether map limits and/or the line limits cross the longitude jump or
% not.
 
  if (e &&                                                              ...
      ((a&&( (b&&c&&d) || (~b&&(c||d)))) || (~a&&(~b||( b&&(c||d)))))   ...
      && g.level<=flaglim)
  
    l=l+1;

    x=C(1:2:end)*1e-6;
    y=C(2:2:end)*1e-6;

% Make things continuous (join edges that cut across 0-meridian)

    dx=diff(x);
    x=x-360*cumsum([x(1)>180; (dx>356) - (dx<-356)]);

% Antarctic is a special case - extend contour to make nice closed polygon
% that doesn't surround the pole.   

    if (g.ver < 15)                      % Range from 0-360
      if (abs(x(1))<1 && abs(y(1)+68.9)<1)
        y=[-89.9;-78.4;y(x<=-180);y(x>-180);    -78.4; -89.9*ones(18,1)];
        x=[  180; 180 ;x(x<=-180)+360;x(x>-180); -180;  [-180:20:160]'];
      end
    else                                 % new version is from -180 to 180
      if (abs(x(1))==180 && y(1)<-76)
        y=[-89.9;-78.4;y;  -78.4; -89.9*ones(19,1)];
        x=[ 180; 180 ;x;   -180;  [-180:20:180]'];
      end
    end  

% First and last point should be the same IF THIS IS A POLYGON
% if the Area=0 then this is a line, and don't add points!
   
    if (g.area > 0)
    
      if (x(end)~=x(1) || y(end)~=y(1))
        x=[x; x(1)];                      % First and last points should
        y=[y; y(1)];                      % be the same
      end
    
% Get correct curve orientation for patch-fill algorithm.
 
      Area2(l)=sum(diff(x).*(y(1:(end-1))+y(2:end))/2); % Area on the page
      Area(l)=g.area/10;                                % Area on the globe

      if (rem(g.level,2) == 0)            % Make lakes (2) and islands (1)
        Area(l)=-abs(Area(l));            % differently oriented
        if (Area2(l) > 0)
          x=x(end:-1:1);
          y=y(end:-1:1);
        end
      else
        if (Area2(l) < 0)
          x=x(end:-1:1);
          y=y(end:-1:1);
        end 
      end
       
    else      % Later on 2 point lines are clipped so we want to avoid that
       
      if (length(x) == 2)
        x=[x(1); mean(x); x(2)];
        y=[y(1); mean(y); y(2)];
      end 
     
    end
 
% Here we try to reduce the number of points and clip them
% close to map limits (if they are too far away this can cause
% problems with some projections as they wrap around into the
% wrong place)
   
    [cx,cy]=clip_to_lims(x,y,m,decfac);

    k(l+1)=k(l)+length(cx)+1;
    ncst(k(l)+1:k(l+1),:)=[cx,cy; NaN NaN];
    type(l)=g.level;
     
% This is a little tricky...the filling algorithm expects data to be in
% the range -180 to 180 deg long. However, there are some land parts that
% cut across this divide so they appear at +190 but not -170. This causes
% problems later... so as a kludge I replicate some of the problematic
% features at 190-360=-170. I cut out the problematic points only, make
% sure the first and last points are the same, and then add this curve to
% ncst.

    if (max(x) > 180)
      ii=find(x>179);ii=[ii;ii(1)];
      [cx,cy]=clip_to_lims(x(ii)-360,y(ii),m,decfac);
      l=l+1;
      k(l+1)=k(l)+length(cx)+1;
      ncst(k(l)+1:k(l+1),:)=[cx,cy;NaN NaN];
    end
    
  end
 
  [cnt,g]=get_gheader(fid);

end

fclose(fid);     

% Get rid of unused part of data matrices.

type(l+1)=type(l);

ncst((k(l+1)+1):end,:)=[];  
Area((l+1):end)=[];
type((l+2):end)=[];
k((l+2):end)=[];    

%--------------------------------------------------------------------------
% Load requested data
%--------------------------------------------------------------------------

C = struct ('lon', [], 'lat', [], 'area', [], 'type', [], 'k', []);

C.lon = ncst(:,1);
C.lat = ncst(:,2);
C.area = Area;
C.type = type;
C.k = k;

return

%--------------------------------------------------------------------------

function [x,y]=clip_to_lims(x,y,m,decfac)

% Look for points outside the lat/lon boundaries, and then decimate them
% by a factor of about 'decfac' (don't get rid of them completely because
% that can sometimes cause problems when polygon edges cross curved map
% edges).
   
tol=.2;   
  
% Do y limits, then x so we can keep corner points.
   
nn=(y>m.tlim+tol) | (y<m.blim-tol);

% Keep one extra point when crossing limits, also the beginning/end point.

nn=logical(nn-min(1,([0;diff(nn)]>0)+([diff(nn);0]<0)));
nn([1 end])=0;

% decimate vigorously

nn=nn & rem(1:length(nn),decfac)'~=0;
x(nn)=[];
y(nn)=[];
         
if (m.rlim > m.llim)                                     % no wraparound
                                % sections of line outside lat/long limits
  nn=(x>m.rlim+tol | x<m.llim-tol) & y<m.tlim & y>m.blim;
else                                                     % wraparound case
  nn=(x>m.rlim+tol & x<m.llim-tol) & y<m.tlim & y>m.blim;
end

nn=logical(nn-min(1,([0;diff(nn)]>0)+([diff(nn);0]<0)));nn([1 end])=0;
nn=nn & rem(1:length(nn),decfac)'~=0;
x(nn)=[];y(nn)=[];
   
% Move all points "near" to map boundaries.
% I'm not sure about the wisdom of this - it might be better to clip
% to the boundaries instead of moving. Hmmm. 
      
y(y>m.tlim+tol)=m.tlim+tol;
y(y<m.blim-tol)=m.blim-tol;

if (m.rlim > m.llim)          % Only clip long bdys if I can tell I'm on
  x(x>m.rlim+tol)=m.rlim+tol; % the right  or left (i.e. not in wraparound
  x(x<m.llim-tol)=m.llim-tol; % case)
end   
   
return

%--------------------------------------------------------------------------

function [cnt,g]=get_gheader(fid)

% Reads the gshhs file header
% 
% A bit of code added because header format changed with version 1.3.
%
% 17/Sep/2008 - added material to handle latest GSHHS version.
%
% For version 1.1 this is the header (9*4 = 36 bytes long)
%
% int id                         Unique polygon id number, starting at 0
% int n                          Number of points in this polygon
% int level                      1 land, 2 lake, 3 island_in_lake,
%                                4 pond_in_island_in_lake
% int west, east, south, north   min/max extent in micro-degrees
% int area                       Area of polygon in 1/10 km^2
% short int greenwich            Greenwich is 1 if Greenwich is crossed
% short int source               0 = CIA WDBII, 1 = WVS
%
% For version 1.3 of GMT format was changed to this ( 10*4 = 40 bytes long)
%
% int id                         Unique polygon id number, starting at 0
% int n                          Number of points in this polygon
% int level                      1 land, 2 lake, 3 island_in_lake,
%                                4 pond_in_island_in_lake
% int west, east, south, north   min/max extent in micro-degrees
% int area                       Area of polygon in 1/10 km^2
% int version                    Polygon version, set to 3
% short int greenwich            Greenwich is 1 if Greenwich is crossed
% short int source               0 = CIA WDBII, 1 = WVS
%
% For version 1.4, we have (8*4 = 32 bytes long)
%
% int id                         Unique polygon id number, starting at 0
% int n                          Number of points in this polygon
% int level                      1 land, 2 lake, 3 island_in_lake,
%                                4 pond_in_island_in_lake
% int flag                       level + version << 8 + greenwich << 16 +
%                                                          source << 24 
% int west, east, south, north   min/max extent in micro-degrees
% int area                       Area of polygon in 1/10 km^2
% int version                    Set to 4 for GSHHS version 1.4
% short int greenwich            Greenwich is 1 if Greenwich is crossed
% short int source               0 = CIA WDBII, 1 = WVS
%
% For version 2.0 it all changed again, we have (11*4 = 44 bytes)
%
% int id                         Unique polygon id number, starting at 0
% int n                          Number of points in this polygon
% int flag                       level + version << 8 + greenwich << 16 +
%                                            source << 24 + river << 25
%
%     flag contains 5 items, as follows:
%  
%     * low byte:      level = flag & 255: Values: 1 land,
%                                                  2 lake,
%                                                  3 island_in_lake,
%                                                  4 pond_in_island_in_lake
%
%                      For Antarctic starting with 2.3.0 ice-front=5 and
%                      grounding line=6
%
%                      For border database: 1=country, 2=state/province
%
%     * 2nd byte:      version = (flag >> 8) & 255:
%                      Values: Should be 7  for GSHHS release 7 (v 2.0)
%                                        12 for v 2.2
%                                        15 for v 2.3.6
%
%     * 3rd byte:      Greenwich = (flag >> 16) & 1:
%                      Values: Greenwich is 1 if Greenwich is crossed
%
%     * 4th byte:      source = (flag >> 24) & 1:
%                      Values: 0 = CIA WDBII, 1 = WVS
%
%     * 5th byte:      river = (flag >> 25) & 1:
%                      Values: 0 = not set, 1 = river-lake and level = 2
%
% int west, east, south, north   min/max extent in micro-degrees
% int area                       Area of polygon in 1/10 km^2 */
% int area_full                  Area of original full-resolution polygon
%                                  in 1/10 km^2
% int container                  Id of container polygon that encloses this
%                                  polygon (-1 if none)
% int ancestor                   Id of ancestor polygon in the full
%                                  resolution set that was the source of
%                                  this polygon (-1 if none) 

% Now, in the calling code I have to use A(2),A(3),A(5-7), A(8), A(9) from
% original.

[A,cnt]=fread(fid,8,'int32');

if (cnt < 8)                               % This gets triggered by the EOF
  g=[];
  return
end
  
g.ver=bitand(bitshift(A(3),-8),255);

if (g.ver == 0)              % then its an old version, fake it to look new

% This works for version 1.2, but not 1.3.

  [A2, cnt2]=fread(fid,1,'int32');
  A=[A; A2];

% We have version 1.3, this would be one of 0, 1, 65535, 65536 in v 1.2
  
  if ((cnt+cnt2)==9 && (A(9)==3))
    A2=fread(fid,1,'int32');                   % Read one more byte
    A(9)=A2;                      % Easiest way not to break existing code
  end                             % one of 0, 1, 65536, 65537

% After v2.0 some more bytes around for original area and container
% and ancerstor IDs.
  
elseif (g.ver >= 7)
  
  A2=fread(fid,3,'int32');

end

% Newest versions.

g.id=A(1);
g.N=A(2);
g.level=bitand(A(3),255);
g.greenwich=bitand(bitshift(A(3),-16),255)*65536;
g.source=bitand(bitshift(A(3),-24),255);
A(3)=g.level;
A(9)=g.greenwich;

g.extentW=A(4);
g.extentE=A(5);
g.extentS=A(6);
g.extentN=A(7);
g.area=A(8);
 
return
