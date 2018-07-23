function [C]=r_gshhs(Llon, Rlon, Blat, Tlat, Fname)

%
% R_GSHHS:  Read and extract coastline data from GSHHS dataset
%
% [Cst,k,C.type,Area]=r_gshhs(Llon, Rlon, Blat, Tlat, Fname)
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
%     Llon      Left   corner longitude (West values are negative).
%     Rlon      Right  corner longitude (West values are negative).
%     Blat      Bottom corner latitude (South values are negative).
%     Tlat      Top    corner latitude (South values are negative).
%     Fname     GSHHS database file name (string).
%
%  On Output:
%
%     C         Read coastline data (structure array):
%                 C.lon  => longitude of closed polygons separated by NaNs.
%                 C.lat  => latitude  of closed polygons separated by NaNs.
%                 C.area => polygon areas.
%                 C.type => polygon type:
%                           C.type=1 => land,
%                           C.type=2 => lake,
%                           C.type=3 => island_in_lake,
%                           C.type=4 => pond_in_island_in_lake
%

% svn $Id: r_gshhs.m 895 2018-02-11 23:15:37Z arango $
%=========================================================================%
%  Copyright (c) 2002-2018 The ROMS/TOMS Group                            %
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

mllim=rem(Llon+360+180,360)-180;
mrlim=rem(Rlon+360+180,360)-180;
mblim=Blat;
mtlim=Tlat;

%--------------------------------------------------------------------------
%  Initialize output data: "decfac" is for decimation of areas outside the
%  lat/lon boundarys.
%--------------------------------------------------------------------------

if (strfind(Fname,'gshhs_f')),         % Full resolution database
  C.lon=NaN+zeros(492283,1);
  C.lat=NaN+zeros(492283,1);
  C.area=zeros(41520,1);
  k=ones(41521,1);
  decfac=12500;
elseif (strfind(Fname,'gshhs_h')),     % High resolution database
  C.lon=NaN+zeros(492283,1);
  C.lat=NaN+zeros(492283,1);
  C.area=zeros(41520,1);
  k=ones(41521,1);
  decfac=2500;
elseif (strfind(Fname,'gshhs_i')),     % Intermediate resolution database
  C.lon=NaN+zeros(492283,1);
  C.lat=NaN+zeros(492283,1);
  C.area=zeros(41520,1);
  k=ones(41521,1);
  decfac=500;
elseif (strfind(Fname,'gshhs_l')),     % Low resolution database
  C.lon=NaN+zeros(101023,1);
  C.lat=NaN+zeros(101023,1);
  C.area=zeros(10768,1);
  k=ones(10769,1);
  decfac=100;
elseif (strfind(Fname,'gshhs_c')),     % Crude resolution databse
  C.lon=NaN+zeros(14872,1);
  C.lat=NaN+zeros(14872,1);
  C.area=zeros(1868,1);
  k=ones(1869,1);
  decfac=20;
else
  error(['READ_GSHHS - unable to process database: ',Fname]);
end

%---------------------------------------------------------------------------
%   Read in GSHHS database.
%---------------------------------------------------------------------------

fid=fopen(Fname,'r','ieee-be');

if (fid==-1),
  error(['read_GSHHS - unable to open file: ',Fname]);
end

Area2=C.area;

%  Read in header (A) for first poligon:
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

[A,cnt]=fread(fid,9,'int32');

l=0;

while cnt > 0,

%  Read all points in the current segment. 

  Cdat=fread(fid,A(2)*2,'int32');
 
  a=rlim > llim;                % Map limits cross longitude jump? (a==1 is no)
  b=A(9) < 65536;               % Cross boundary? (b==1 if no).
  c=llim < rem(A(5)+360e6,360e6); 
  d=rlim > rem(A(4)+360e6,360e6);
  e=tlim > A(6) & blim<A(7);
 
%  This test checks whether the lat/long box containing the line overlaps that
%  of the map. There are various cases to consider, depending on whether map
%  limits and/or the line limits cross the longitude jump or not.
 
  if (e && (a && ((b&&c&&d) || (~b&&(c||d))) || (~a && (~b || (b&&(c||d)))))),
 
    l=l+1;

    x=Cdat(1:2:end)*1e-6;
    y=Cdat(2:2:end)*1e-6;

%  Make things continuous (join edges that cut across 0-meridian)

    dx=diff(x);
    x=x-360*cumsum([x(1)>180; (dx>356) - (dx<-356)]);

%  Antarctic is a special case - extend contour to make nice closed polygon
%  that doesn't surround the pole.   

    if ( abs(x(1))<1 && abs(y(1)+68.9)<1 ),
      y=[-89.9; -78.4; y(x<=-180);     y(x>-180); -78.4; -89.9*ones(18,1)];
      x=[  180;   180; x(x<=-180)+360; x(x>-180);  -180; [-180:20:160]'];
    end

%  First and last point should be the same.
   
    if ( x(end)~=x(1) || y(end)~=y(1) ),
      x=[x; x(1)];
      y=[y; y(1)];
    end

%  Get correct curve orientation for patch-fill algorithm.
   
    Area2(l)=sum( diff(x).*(y(1:(end-1))+y(2:end))/2 );
    C.area(l)=A(8)/10;

    if (rem(A(3),2) == 0),
      C.area(l)=-abs(C.area(l)); 
      if (Area2(l) > 0),
        x=x(end:-1:1);
        y=y(end:-1:1);
      end,
    else
      if (Area2(l) < 0),
        x=x(end:-1:1);
        y=y(end:-1:1);
      end, 
    end,

%  Here we try to reduce the number of points.
   
    xflag=0;
    if (max(x) > 180),  % First, save original curve for later if we anticipate
      sx=x;             % a 180-problem.
      sy=y;         
      xflag=1;
    end;
   
%  Look for points outside the lat/long boundaries, and then decimate them by
%  a factor of about 'decfac' (don't get rid of them completely because that
%  can sometimes cause problems when polygon edges cross curved map edges).
   
    tol=0.2;   
  
%  Do "y" limits, then "x" so we can keep corner points.
   
    nn=(y > mtlim+tol) | (y < mblim-tol);

%  Keep one extra point when crossing limits, also the beginning/end point.

    nn=logical(nn-([0;diff(nn)]>0)-([diff(nn);0]<0));
    nn([1 end])=0;

%  Decimate vigorously.

    nn=nn & rem(1:length(nn),decfac)'~=0;
    x(nn)=[];
    y(nn)=[];
         
    if (mrlim > mllim),                                   % no wraparound case
      nn=(x>mrlim+tol | x<mllim-tol) & y<mtlim & y>mblim;
    else                                                  % wraparound case
      nn=(x>mrlim+tol & x<mllim-tol ) & y<mtlim & y>mblim;
    end;
    nn=logical(nn-([0;diff(nn)]>0)-([diff(nn);0]<0));
    nn([1 end])=0;
    nn=nn & rem(1:length(nn),decfac)'~=0;
    x(nn)=[];
    y(nn)=[];
   
%  Move all points "near" to map boundaries.  I'm not sure about the wisdom
%  of this - it might be better to clip to the boundaries instead of moving.
      
    y(y>mtlim+tol)=mtlim+tol;
    y(y<mblim-tol)=mblim-tol;

%  Only clip longitude boundarys if I can tell I'm on the right or left
%  (that is, not in wraparound case).
 
    if (mrlim > mllim),
      x(x > mrlim+tol)=mrlim+tol;
      x(x < mllim-tol)=mllim-tol;
    end;   
   
    k(l+1)=k(l)+length(x)+1;
    C.lon(k(l)+1:k(l+1)-1)=x;  C.lon(k(l+1))=NaN;
    C.lat(k(l)+1:k(l+1)-1)=y;  C.lat(k(l+1))=NaN;
    C.type(l)=A(3);
   
%  This is a little tricky...the filling algorithm expects data to be in the
%  range -180 to 180 deg long. However, there are some land parts that cut
%  across this divide so they appear at +190 but not -170. This causes
%  problems later... so as a kludge I replicate some of the problematic
%  features at 190-360=-170. Small islands are just duplicated, for the
%  Eurasian landmass I just clip off the eastern part.
   
    if (xflag),
      l=l+1;
      C.area(l)=C.area(l-1); 
      if (abs(C.area(l)) > 1e5),
        nn=find(sx>180);
        nn=[nn; nn(1)];
        k(l+1)=k(l)+length(nn)+1;
        C.lon(k(l)+1:k(l+1)-1)=sx(nn)-360;
        C.lat(k(l)+1:k(l+1)-1)=sy(nn);

      else,                              % repeat the island at the other edge.
        k(l+1)=k(l)+length(sx)+1;
        C.lon(k(l)+1:k(l+1)-1)=sx-360;
        C.lat(k(l)+1:k(l+1)-1)=sy;
      end;
      C.lon(k(l+1))=NaN;
      C.lat(k(l+1))=NaN;
      C.type(l)=A(3);
    end;

  end;

%  Read header for next polygon.
 
  [A,cnt]=fread(fid,9,'int32');

end,

fclose(fid);

C.type=C.type';
C.type(l+1)=C.type(l);

%  Get rid of unused part of data matrices.

C.lon((k(l+1)+1):end)=[];
C.lat((k(l+1)+1):end)=[];
C.area((l+1):end)=[];
C.type((l+2):end)=[];

return
