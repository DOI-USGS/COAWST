function [Istr,Iend,Jstr,Jend] = sample_grid(XD,YD,XR,YR,varargin)

%
% SAMPLE_GRID:  Gets donor grid tight indices containing receiver grid
%
% [Istr,Iend,Jstr,Jend] = sample_grid(XD,YD,XR,YR,offset,plt);
%
% This function computes the donor grid indices range of the polygon
% that tightly contains the receiver grid.  This is done to sample the
% donor grid to accelerate the interpolation of fields for the receiver
% grid.
%
% On Input:
%
%    XD            Donor grid X-coordinates (2D array)
%    YD            Donor grid Y-coordinates (2D array)
%    XR            Receiver grid X-coordinates (2D array)
%    YR            Receiver grid Y-coordinates (2D array)
%    offset        Number of extra points used to sample the
%                    donor grid so it is large enough to contain
%                    the receiver grid  (default 5)
%    plt           Switch to plot donor and receiver grids
%                    (default false)
%
% On Output:
%
%    Istr          Donor grid starting I-index for sampling
%    Iend          Donor grid ending   I-index for sampling
%    Jstr          Donor grid starting J-index for sampling
%    Jend          Donor grid ending   J-index for sampling
%

% svn $Id: sample_grid.m 738 2014-10-14 21:49:14Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set optional arguments.

offset = 5;
plt = false;

switch numel(varargin)
  case 1
    offset = varargin{1};
  case 2
    offset = varargin{1};
    plt = varargin{2};
end

%--------------------------------------------------------------------------
%  Extract donor grid polygon.
%--------------------------------------------------------------------------

[ImD,JmD] = size(XD);
IstrD = 1;
IendD = ImD;
JstrD = 1;
JendD = JmD;

XboxD = [squeeze(XD(IstrD:IendD,JstrD));                              ...
         squeeze(XD(IendD,JstrD+1:JendD))';                           ...
         squeeze(flipud(XD(IstrD:IendD-1,JendD)));                    ...
         squeeze(fliplr(XD(IstrD,JstrD:JendD-1)))'];

YboxD = [squeeze(YD(IstrD:IendD,JstrD));                              ...
         squeeze(YD(IendD,JstrD+1:JendD))';                           ...
         squeeze(flipud(YD(IstrD:IendD-1,JendD)));                    ...
         squeeze(fliplr(YD(IstrD,JstrD:JendD-1)))'];

%--------------------------------------------------------------------------
%  Extract receiver grid polygon.
%--------------------------------------------------------------------------

[ImR,JmR] = size(XR);
IstrR = 1;
IendR = ImR;
JstrR = 1;
JendR = JmR;

XboxR = [squeeze(XR(IstrR:IendR,JstrR));                                ...
         squeeze(XR(IendR,JstrR+1:JendR))';                             ...
         squeeze(flipud(XR(IstrR:IendR-1,JendR)));                      ...
         squeeze(fliplr(XR(IstrR,JstrR:JendR-1)))'];

YboxR = [squeeze(YR(IstrR:IendR,JstrR));                                ...
         squeeze(YR(IendR,JstrR+1:JendR))';                             ...
         squeeze(flipud(YR(IstrR:IendR-1,JendR)));                      ...
         squeeze(fliplr(YR(IstrR,JstrR:JendR-1)))'];

%--------------------------------------------------------------------------
%  Find donor indices containing the receiver grid.
%--------------------------------------------------------------------------

[in,on] = inpolygon(XD,YD,XboxR,YboxR);

%--------------------------------------------------------------------------
%  Set donor grid sampling indices.
%--------------------------------------------------------------------------

[J,I] = meshgrid(1:1:JmD, 1:1:ImD);

I(~in) = NaN;
J(~in) = NaN;

Istr = min(I(:))-offset;
if (isnan(Istr) || Istr < 1),
  Istr = 1;
end

Iend = max(I(:))+offset;
if (isnan(Iend) || Iend > ImD),
  Iend = ImD;
end

Jstr = min(J(:))-offset;
if (isnan(Jstr) || Jstr < 1),
  Jstr = 1;
end

Jend = max(J(:))+offset;
if (isnan(Jend) || Jend > JmD),
  Jend = JmD;
end

%--------------------------------------------------------------------------
%  If requested, plot donor and receiver grids.
%--------------------------------------------------------------------------

if (plt),

  X = XD(Istr:1:Iend,Jstr:1:Jend);
  Y = YD(Istr:1:Iend,Jstr:1:Jend);

  figure;
  plot(XboxD,YboxD,'r.',XboxR,YboxR,'b.');
  title(['Original Donor and Receiver Grids']);
  
  figure;
  pcolor(X,Y,ones(size(X)));
  hold on;
  plot(XboxR,YboxR,'b.');
  title(['Sampled Donor and Receiver Grids']);
  hold off

end

return
