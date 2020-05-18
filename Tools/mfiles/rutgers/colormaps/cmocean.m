function cmap = cmocean(ColormapName,varargin) 

% CMOCEAN returns perceptually-uniform colormaps created by Kristen Thyng. 
% 
%% Syntax 
% 
%  cmocean 
%  cmap = cmocean('ColormapName') 
%  cmap = cmocean('-ColormapName') 
%  cmap = cmocean(...,NLevels)
%  cmap = cmocean(...,'pivot',PivotValue) 
%  cmap = cmocean(...,'negative') 
%  cmocean(...)
% 
%% Description 
% 
% cmocean without any inputs displays the options for colormaps. 
% 
% cmap = cmocean('ColormapName') returns a 256x3 colormap. ColormapName can be any of 
% of the following: 
% 
%          SEQUENTIAL:                DIVERGING: 
%          'thermal'                  'balance'
%          'haline'                   'delta'
%          'solar'                    'diff'
%          'ice'                      'curl'
%          'gray'                     'tarn'
%          'oxy'
%          'deep'                     CONSTANT LIGHTNESS:
%          'dense'                    'phase'
%          'algae'
%          'matter'
%          'turbid'
%          'speed'
%          'amp'
%          'tempo'
%          'rain'
%
%          OTHERS:
%          'cividis'		      
%          'inferno'  
%          'magma'
%          'parula'    (Matlab default)
%          'plasma'
%          'topo'
%          'viridis'
%
% cmap = cmocean('-ColormapName') a minus sign preceeding any ColormapName flips the
% order of the colormap. 
%
% cmap = cmocean(...,NLevels) specifies a number of levels in the colormap.  Default
% value is 256. 
%
% cmap = cmocean(...,'pivot',PivotValue) centers a diverging colormap such that white 
% corresponds to a given value and maximum extents are set using current caxis limits. 
% If no PivotValue is set, 0 is assumed. Early versions of this function used 'zero'
% as the syntax for 'pivot',0 and the old syntax is still supported. 
%
% cmap = cmocean(...,'negative') inverts the lightness profile of the colormap. This can be 
% useful particularly for divergent colormaps if the default white point of divergence
% gets lost in a white background. 
% 
% cmocean(...) without any outputs sets the current colormap to the current axes.  
% 
%% Examples:
%
% Using this sample plot: 
% 
%   imagesc(peaks(1000)+1)
%   colorbar
% 
% Set the colormap to 'algae': 
% 
%   cmocean('algae') 
% 
% Same as above, but with an inverted algae colormap: 
% 
%   cmocean('-algae')
% 
% Set the colormap to a 12-level 'solar': 
% 
%   cmocean('solar',12)
% 
% Get the RGB values of a 5-level thermal colormap: 
% 
%   RGB = cmocean('thermal',5)
% 
% Some of those values are below zero and others are above. If this dataset represents
% anomalies, perhaps a diverging colormap is more appropriate: 
% 
%   cmocean('balance') 
% 
% It's unlikely that 1.7776 is an interesting value about which the data values 
% diverge.  If you want to center the colormap on zero using the current color 
% axis limits, simply include the 'pivot' option:  
% 
%   cmocean('balance','pivot',0) 
%
%% Author Info:
%
% This function was written by Chad A. Greene of the Institute for Geophysics at the 
% University of Texas at Austin (UTIG), June 2016, using colormaps created by Kristen
% Thyng of Texas A&M University, Department of Oceanography. More information on the
% cmocean project can be found at http://matplotlib.org/cmocean/. 
% 
%% Citing this colormap: 
%
% If you find an occasion to cite these colormaps for any reason, or if you just want
% some nice beach reading, check out the following paper from the journal Oceanography: 
% 
% Kristen M. Thyng, Chad A. Greene, Robert D. Hetland, Heather M. Zimmerle, and Steven
% F. DiMarco. True colors of oceanography: Guidelines for effective and accurate colormap
% selection. Oceanography, September 2016, http://dx.doi.org/10.5670/oceanog.2016.66
% 
% See also colormap,  caxis,    colorspace,  
%          cm_algae,  cm_amp,   cm_balance, cm_curl,  cm_deep,    cm_delta,
%          cm_dense,  cm_gray,  cm_haline,  cm_ice,   cm_matter,  cm_oxy
%          cm_phase,  cm_solar, cm_speed,   cm_tempo, cm_thermal, cm_turbid,
%          inferno,   magma,    parula,     plasma,   cividis,    viridis
  
% svn $Id: cmocean.m 996 2020-01-10 04:28:56Z arango $

% Display colormap options: 

if nargin==0
   figure('menubar','none','numbertitle','off','Name', ...
          'cmocean options:','pos',[10 10 364 751])
   axes('pos',[0 0 1 1])
   image(imread('cmocean.png')); 
   axis image off
   return
end

% Error checks: 

assert(isnumeric(ColormapName)==0,'Input error: ColormapName must be a string.')

% Set defaults: 

NLevels = 256; 
autopivot = false; 
PivotValue = 0; 
InvertedColormap = false; 

% Parse inputs: 

% Does user want to flip the colormap direction? 

dash = regexp(ColormapName,'-'); 
if any(dash) 
   InvertedColormap = true; 
   ColormapName(dash) = []; 
end

% Forgive the British: 

if strncmpi(ColormapName,'grey',4); 
   ColormapName = 'gray'; 
end

% Does the user want a "negative" version of the colormap
% (with an inverted lightness profile)? 

tmp = strncmpi(varargin,'negative',3); 
if any(tmp) 
   negativeColormap = true; 
   varargin = varargin(~tmp); 
else
   negativeColormap = false; 
end

% Does the user want to center a diverging colormap on a specific value? 
% This parsing support original 'zero' syntax and current 'pivot' syntax. 

tmp = any([strncmpi(varargin,'pivot',3) strncmpi(varargin,'zero',3)]); 
if any(tmp) 
   autopivot = true; 
   try
      if isscalar(varargin{find(tmp)+1})
         PivotValue = varargin{find(tmp)+1}; 
         tmp(find(tmp)+1) = 1; 
      end
   end
   varargin = varargin(~tmp); 
end

% Has user requested a specific number of levels? 

tmp = isscalar(varargin); 
if any(tmp) 
   NLevels = varargin{tmp}; 
end

% Load RGB values and interpolate to NLevels: 

switch lower(ColormapName(1:3))
   case 'alg'
      cmap = cm_algae;
   case 'amp' 
      cmap = cm_amp;
   case 'bal' 
      cmap = cm_balance;
   case 'civ' 
      cmap = cividis;
   case 'cur'
      cmap = cm_curl;
   case 'dee' 
      cmap = cm_deep;
   case 'del'
      cmap = cm_delta;
   case 'den' 
      cmap = cm_dense;
   case 'dif' 
      cmap = cm_diff;
   case 'gra'
      cmap = cm_gray;
   case 'hal' 
      cmap = cm_haline;
   case 'ice'
      cmap = cm_ice;
   case 'inf' 
      cmap = inferno;
   case 'mag' 
      cmap = magma;
   case 'mat' 
      cmap = cm_matter;
   case 'oxy' 
      cmap = cm_oxy;
   case 'par' 
      cmap = parula;          % Matlab default
   case 'pha' 
      cmap = cm_phase;
   case 'pla' 
      cmap = plasma;
   case 'rai' 
      cmap = cm_rain;
   case 'sol' 
      cmap = cm_solar;
   case 'spe' 
      cmap = cm_speed;
   case 'tar' 
      cmap = cm_tarn;
   case 'tem' 
      cmap = cm_tempo;
   case 'the' 
      cmap = cm_thermal;
   case 'top' 
      cmap = cm_topo;
   case 'tur' 
      cmap = cm_turbid;
   case 'vir' 
      cmap = viridis;
   otherwise 
      error('Unrecognized colormap name.') 
end

%  Transform color image representation:

if negativeColormap
   
   % Convert RGB to LAB colorspace: 

   LAB = colorspace('RGB->LAB',cmap); 

   % Operate on the lightness profile: 

   L = LAB(:,1); 

   % Flip the lightness profile and set the lowest point to black:

   L = max(L) - L; 

   % Stretch the lightness profile to make the lightest bits 95% white.
   % (Going 100% white would make the ends of a divergent profile
   % impossible to distinguish.)
   
   L = L*(95/max(L)); 

   % Make a new LAB matrix: 

   LAB = [L LAB(:,2:3)]; 
   
   % Convert LAB back to RGB: 
   
   cmap = colorspace('LAB->RGB',LAB); 
end

% Interpolate if necessary: 

if NLevels~=size(cmap,1); 
   R = interp1((1:size(cmap,1)),cmap(:,1),linspace(1,size(cmap,1),NLevels)'); 
   G = interp1((1:size(cmap,1)),cmap(:,2),linspace(1,size(cmap,1),NLevels)'); 
   B = interp1((1:size(cmap,1)),cmap(:,3),linspace(1,size(cmap,1),NLevels)'); 
   cmap = [R G B]; 
end

% Invert the colormap if requested by user: 

if InvertedColormap
   cmap = flipud(cmap); 
end

%% Adjust values to current caxis limits? 

if autopivot
   clim = caxis; 
   maxval = max(abs(clim-PivotValue)); 
   
   R = interp1( linspace(-maxval,maxval,size(cmap,1))'+PivotValue ,cmap(:,1), ...
                linspace(clim(1),clim(2),size(cmap,1))'); 
   G = interp1( linspace(-maxval,maxval,size(cmap,1))'+PivotValue ,cmap(:,2), ...
                linspace(clim(1),clim(2),size(cmap,1))'); 
   B = interp1( linspace(-maxval,maxval,size(cmap,1))'+PivotValue ,cmap(:,3), ...
                linspace(clim(1),clim(2),size(cmap,1))'); 
   cmap = [R G B]; 
end

% Clean up

if nargout==0
   colormap(gca,cmap) 
   clear cmap  
end

return


