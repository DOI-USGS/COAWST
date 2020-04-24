function [C]=x_gshhs(Llon, Rlon, Blat, Tlat, C, cliptype)

%
% X_GSHHS:  Process extracted GSHHS coastline data
%
% [C]=x_gshhs(Llon, Rlon, Blat, Tlat, C, cliptype)
%
%  This function process coastline data read from GSHHS database.
%
%  Adapted from Rich Pawlowicz "mu_coast" function.
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
%     C         Read coastline data (structure array):
%                 C.lon  => longitude of closed polygons separated by NaNs
%                 C.lat  => latitude  of closed polygons separated by NaNsa
%                 C.area => polygon areas.
%                 C.type => polygon type:
%                           C.type=1 => land,
%                           C.type=2 => lake,
%                           C.type=3 => island_in_lake,
%                           C.type=4 => pond_in_island_in_lake
%                 C.k    => number of points in polygon
%
%     cliptype  Clip type:
%                 'on'    - replaces points outside with NaN, but
%                           interpolates to provide points right on the
%                           border.
%                 'patch' - replaces points outside with nearest border
%                           point.
%                 'point' - does no interpolation (just checks in/out)
%
%  On Output:
%
%     C         Updated coastline data (structure array).
%

% svn $Id: x_gshhs.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  Handle wrap-arounds.
%--------------------------------------------------------------------------

if (Rlon < -180)
  C.lon=C.lon-360;
elseif (Llon > 180)
  C.lon=C.lon+360;
elseif (Llon < -180)
  C.area=[C.area; C.area];
  C.lon=[C.lon; C.lon(2:end)-360];
  C.lat=[C.lat; C.lat(2:end)];

% This is kinda kludgey - but sometimes adding all these extra points
% causes wrap-around in the conic projection, so we want to limit the
% longitudes to the range needed. However, we don't just clip them to
% min long because that can cause problems in trying to decide which way
% curves are oriented when doing the fill algorithm. So instead I sort
% of crunch the scale, preserving topology.

  nn=C.lon < Llon;
  C.lon(nn)=(C.lon(nn)-Llon)/100 + Llon;
elseif (Rlon > 180)
  C.area=[C.area; C.area];
  C.lon=[C.lon; C.lon(2:end)+360];
  C.lat=[C.lat; C.lat(2:end)];

  nn=C.lon > Rlon;
  C.lon(nn)=(C.lon(nn)-Rlon)/100 + Rlon;
end

%--------------------------------------------------------------------------
%  Clip out-of-range values.
%--------------------------------------------------------------------------

lon=C.lon;
lat=C.lat;

ind=find(isnan(lon));
type=-2.*ones(size(lon));

[lon,lat,type]=m_clip(cliptype,lon,Llon,lon<Llon,lat,type);
[lon,lat,type]=m_clip(cliptype,lon,Rlon,lon>Rlon,lat,type);
[lat,lon,type]=m_clip(cliptype,lat,Blat,lat<Blat,lon,type);
[lat,lon,type]=m_clip(cliptype,lat,Tlat,lat>Tlat,lon,type);

%--------------------------------------------------------------------------
%  Remove redundant NaNs.
%--------------------------------------------------------------------------

switch cliptype
  case 'on'
    ind1=find(isnan(lon));
    ic=0;
    for i=2:length(ind1)
      if ((ind1(i)-ind1(i-1)) == 1)
        ic=ic+1;
        ind2(ic)=ind1(i);
      end
    end
    if (~isempty(ind2))
      ind2=ind2';
      lon(ind2)=[];
      lat(ind2)=[];
      type(ind2)=[];
    end
    ind=find(type == -2);
    type(ind)=[];
    C.type=type;
end

C.lon=lon;
C.lat=lat;

return

function [Xc,Yc,type]=m_clip(cliptype, X, Xedge, index, Y, type)

%
% M_CLIP:  Clip extrated coastline data
%
% [Xc,Yc]=m_clip(cliptype, X, Xedge, index, Y)
%
%  This performs clipping of data. Columns of points are assumed to be
%  lines; the first points outside the clip area are recomputed to lie
%  on the edge of the region; others are converted to either NaN or edge
%  edge points depending on CLIPTYPE. In general m_clip will be called
%  four times, since it solves a 1-edge problem.
%
%  Adapted from Rich Pawlowicz "m_clip" function.
%
%  On Input:
%
%     cliptype  Clip type:
%                 'on'    - replaces points outside with NaN, but
%                           interpolates to provide points right on the
%                           border.
%                 'patch' - replaces points outside with nearest border
%                           point.
%                 'point' - does no interpolation (just checks in/out)
%     index     0 for points inside the clip region, 1 for those outside.
%     type      polygon type.
%

Xc=X;
Yc=Y;

if (~strcmp(cliptype,'point'))

%  Find regions where we suddenly come into the area
%  (indices go from 1 to 0)

  [i,j]=find(diff(index)==-1);

  if (any(i))

    I=i+(j-1)*size(X,1);                 % 1-d addressing

% Linearly interpolate to boundary

    bt=(X(I+1)-X(I));
    ibt=abs(bt)<5*eps;

% In these cases the delta(Y) also=0, so we just want to avoid /0 warnings.

    if (any(ibt))
      bt(ibt)=1*eps;
    end 

    Yc(I)=Y(I)+(Xedge-X(I)).*(Y(I+1)-Y(I))./bt;
    Yc(I(ibt))=(Y(I(ibt))+Y(I(ibt)+1))/2;
    Xc(I)=Xedge;
    index(I(isfinite(Yc(I))))=0;

  end

%  Find regions where we suddenly come out of the area
%  (indices go from 0 to 1)

  [i,j]=find(diff(index)==1);

  if any(i)

    I=i+(j-1)*size(X,1);

    bt=(X(I+1)-X(I));
    ibt=abs(bt)<5*eps;

% In these cases the delta(Y) also=0, so we just want to avoid /0 warnings.

    if (any(ibt))
      bt(ibt)=eps;
    end
  
    Yc(I+1)=Y(I)+(Xedge-X(I)).*(Y(I+1)-Y(I))./bt;
    Yc(I(ibt)+1)=(Y(I(ibt))+Y(I(ibt)+1))/2;
    Xc(I+1)=Xedge;
    index(I(isfinite(Yc(I+1)))+1)=0;
  end

end

switch cliptype
  case {'on','point'}
    Xc(index)=NaN;
    Yc(index)=NaN;
    type(index)=-1;
  case 'patch'
    Xc(index)=Xedge;
end

return
