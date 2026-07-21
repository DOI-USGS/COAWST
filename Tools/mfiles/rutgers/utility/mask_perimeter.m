function [X,Y] = mask_perimeter (Mask, varargin)

% MASK_PERIMETER:  Finds the enclosed polygon perimeter of a mask.
%
% [X,Y] = mask_perimeter(Mask, MaskValue, Lplot)
%
% Given a 2D mask, this function fine the enclosed polygom perimeter for
% all associated points.
%
% On Input:
%
%    Mask        2D mask array (numbers or logicals)
%
%    MaskValue   Mask value to process (Optional)
%                  (default = 1)
%
%    Lplot       Switch to plot mask and perimeters (Optional)
%                  (default = false)
%
% On Output:
%
%    X           I-index location of mask perimeter (vector)
%
%    Y           J-index location of mask perimeter (vector)
%
% Calls:  
%
%    contourc    Computes contour matrix
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%  

% Initialize.

switch numel(varargin)
  case 0
    MaskValue = 1;
    Lplot = false;
  case 1
    MaskValue = varargin{1};
    Lplot = false;
  case 2
    MaskValue = varargin{1};
    Lplot = varargin{2};
end

% Compute countour matrix.

V = [0.5*MaskValue 0.5*MaskValue];
C = contourc(Mask, V);

% Add NaN in between the segments to separate distinct polygons defining
% each unique land mass.

[x,y] = C2xyz(C);

for k=1:length(x)
  x{k}=[x{k} NaN];
  y{k}=[y{k} NaN];
end

X=cell2mat(y);
Y=cell2mat(x);

% Remove outliers.

[Im,Jm] = size(Mask);

ii = find(X < 1 | X > Im+2);
if ~isempty(ii)
  X(ii) = [];
  Y(ii) = [];
end

jj = find(Y < 1 | Y > Jm+2);
if ~isempty(jj)
  X(jj) = [];
  Y(jj) = [];
end

% Plot mask outline boundaries.

if (Lplot)
  [I,J] = find(Mask == MaskValue);
  
  figure;
  plot(I, J, 'r.', X, Y, 'ko-')
  axis([1 Im 1 Jm]);
  grid on;
  title(['Enclosed Mask Perimeter, MaskValue = ',num2str(MaskValue)]);  
end

return

