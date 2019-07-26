%$$$
%$$$ #undef __PR
%$$$ #include "VARIANT.h"

function c = sw_satAr(S,T)

% SW_SATAr   Satuaration of Ar in sea water
%=========================================================================
% sw_satAr $Id: sw_satAr.m 330 2009-03-10 05:57:42Z arango $
%          Copyright (C) CSIRO, Phil Morgan 1998.
%
% USAGE:  satAr = sw_satAr(S,T)
%
% DESCRIPTION:
%    Solubility (satuaration) of Argon (Ar) in sea water
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (ITS-90)]
%
% OUTPUT:
%   satAr = solubility of Ar  [ml/l]
%
% AUTHOR:  Phil Morgan 97-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)
%
%$$$ #include "disclaimer_in_code.inc"
%
% REFERENCES:
%    Weiss, R. F. 1970
%    "The solubility of nitrogen, oxygen and argon in water and seawater."
%    Deap-Sea Research., 1970, Vol 17, pp721-735.
%

% Modifications
% 99-06-25. Lindsay Pender, Fixed transpose of row vectors.
% 03-12-12. Lindsay Pender, Converted to ITS-90.

% CALLER: general purpose
% CALLEE:

%$$$ #ifdef VARIANT_PRIVATE
%$$$ %***********************************************************
%$$$ %
%$$$ %$Log: sw_satAr.m,v $
%$$$ %Revision 1.1  2003/12/12 04:23:22  pen078
%$$$ %*** empty log message ***
%$$$ %

%$$$ %
%$$$ %***********************************************************
%$$$ #endif

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=2
   error('Must pass 2 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);


% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
   error('S & T must have same dimensions')
end %if

%------
% BEGIN
%------

% convert T to Kelvin
T = 273.15 + T * 1.00024;

% constants for Eqn (4) of Weiss 1970
a1 = -173.5146;
a2 =  245.4510;
a3 =  141.8222;
a4 =  -21.8020;
b1 =   -0.034474;
b2 =    0.014934;
b3 =   -0.0017729;

% Eqn (4) of Weiss 1970
lnC = a1 + a2.*(100./T) + a3.*log(T./100) + a4.*(T./100) + ...
      S.*( b1 + b2.*(T./100) + b3.*((T./100).^2) );

c = exp(lnC);

return


