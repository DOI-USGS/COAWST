function c3515 = sw_c3515()

% SW_C3515   Conductivity at (35,15,0)
%=========================================================================
% SW_c3515  $Id: sw_c3515.m 330 2009-03-10 05:57:42Z arango $
%       %   Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  c3515 = sw_c3515
%
% DESCRIPTION:
%   Returns conductivity at S=35 psu , T=15 C [ITPS 68] and P=0 db).
%
% INPUT: (none)
%
% OUTPUT:
%   c3515  = Conductivity   [mmho/cm == mS/cm]
%
% AUTHOR:  Phil Morgan 93-04-17  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    R.C. Millard and K. Yang 1992.
%    "CTD Calibration and Processing Methods used by Woods Hole
%     Oceanographic Institution"  Draft April 14, 1992
%    (Personal communication)
%

% CALLER: none
% CALLEE: none
%

c3515 = 42.914;

return
