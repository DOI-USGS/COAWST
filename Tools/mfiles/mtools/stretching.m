function [s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid,  ...
                          report)

%
% STRETCHING:  Compute ROMS vertical coordinate stretching function
%
% [s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report)
%
% Given vertical terrain-following vertical stretching parameters, this
% routine computes the vertical stretching function used in ROMS vertical
% coordinate transformation. Check the following link for details:
%
%    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
%
% On Input:
%
%    Vstretching   Vertical stretching function:
%                    Vstretching = 1,  original (Song and Haidvogel, 1994)
%                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
%                    Vstretching = 3,  R. Geyer BBL refinement
%                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)
%    theta_s       S-coordinate surface control parameter (scalar)
%    theta_b       S-coordinate bottom control parameter (scalar)
%    hc            Width (m) of surface or bottom boundary layer in which
%                    higher vertical resolution is required during
%                    stretching (scalar)
%    N             Number of vertical levels (scalar)
%    kgrid         Depth grid type logical switch:
%                    kgrid = 0,        function at vertical RHO-points
%                    kgrid = 1,        function at vertical W-points
%    report        Flag to report detailed information (OPTIONAL):
%                    report = false,   do not report
%                    report = true,    report information
%
% On Output:
%
%    s             S-coordinate independent variable, [-1 <= s <= 0] at
%                    vertical RHO- or W-points (vector)
%    C             Nondimensional, monotonic, vertical stretching function,
%                    C(s), 1D array, [-1 <= C(s) <= 0]
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

s=[];
C=[];

%--------------------------------------------------------------------------
%  Set several parameters.
%--------------------------------------------------------------------------

if (nargin < 6),
  disp(' ');
  disp('*** Error:  STRETCHING - too few arguments.');
  disp(['                     number of supplied arguments: ',          ...
        num2str(nargin)]);
  disp('                     number of required arguments: 6');
  disp(' ');
  return
end,

if (Vstretching < 1 || Vstretching > 4),
  disp(' ');
  disp(['*** Error:  STRETCHING - Illegal parameter Vstretching = '     ...
	num2str(Vstretching)]); 
  disp(' ');
  return
end

if (nargin < 7),
  report=false;
end

Np=N+1;

%--------------------------------------------------------------------------
% Compute ROMS S-coordinates vertical stretching function
%--------------------------------------------------------------------------

% Original vertical stretching function (Song and Haidvogel, 1994).

if (Vstretching == 1),

  ds=1.0/N;
  if (kgrid == 1),
    Nlev=Np;

    lev=0:N;
    s=(lev-N).*ds;
  else
    Nlev=N;
    lev=(1:N)-0.5;
    s=(lev-N).*ds;
  end
  if (theta_s > 0),
    Ptheta=sinh(theta_s.*s)./sinh(theta_s);
    Rtheta=tanh(theta_s.*(s+0.5))./(2.0*tanh(0.5*theta_s))-0.5;
    C=(1.0-theta_b).*Ptheta+theta_b.*Rtheta;
  else
    C=s;
  end

% A. Shchepetkin (UCLA-ROMS, 2005) vertical stretching function.

elseif (Vstretching == 2),

  alfa=1.0;
  beta=1.0;
  ds=1.0/N;
  if (kgrid == 1),
    Nlev=Np;
    lev=0:N;
    s=(lev-N).*ds;
  else
    Nlev=N;
    lev=(1:N)-0.5;
    s=(lev-N).*ds;
  end
  if (theta_s > 0),
    Csur=(1.0-cosh(theta_s.*s))/(cosh(theta_s)-1.0);
    if (theta_b > 0),
      Cbot=-1.0+sinh(theta_b*(s+1.0))/sinh(theta_b);
      weigth=(s+1.0).^alfa.*(1.0+(alfa/beta).*(1.0-(s+1.0).^beta));
      C=weigth.*Csur+(1.0-weigth).*Cbot;
    else
      C=Csur;
    end
  else
    C=s;
  end

%  R. Geyer BBL vertical stretching function.

elseif (Vstretching == 3),

  ds=1.0/N;
  if (kgrid == 1),
    Nlev=Np;
    lev=0:N;
    s=(lev-N).*ds;
  else
    Nlev=N;
    lev=(1:N)-0.5;
    s=(lev-N).*ds;
  end
  if (theta_s > 0),
     exp_s=theta_s;      %  surface stretching exponent
     exp_b=theta_b;      %  bottom  stretching exponent
     alpha=3;            %  scale factor for all hyperbolic functions
    Cbot=log(cosh(alpha*(s+1).^exp_b))/log(cosh(alpha))-1;
    Csur=-log(cosh(alpha*abs(s).^exp_s))/log(cosh(alpha));
    weight=(1-tanh( alpha*(s+.5)))/2;
    C=weight.*Cbot+(1-weight).*Csur;
  else
    C=s;
  end

% A. Shchepetkin (UCLA-ROMS, 2010) double vertical stretching function
% with bottom refinement

elseif (Vstretching == 4),

  ds=1.0/N;
  if (kgrid == 1),
    Nlev=Np;
    lev=0:N;
    s=(lev-N).*ds;
  else
    Nlev=N;
    lev=(1:N)-0.5;
    s=(lev-N).*ds;
  end
  if (theta_s > 0),
    Csur=(1.0-cosh(theta_s.*s))/(cosh(theta_s)-1.0);
  else
    Csur=-s.^2;
  end
  if (theta_b > 0),
    Cbot=(exp(theta_b.*Csur)-1.0)/(1.0-exp(-theta_b));
    C=Cbot;
  else
    C=Csur;
  end

end

% Report S-coordinate parameters.

if (report),
  disp(' ');
  if (Vstretching == 1),
    disp(['Vstretching = ',num2str(Vstretching),                        ...
          '   Song and Haidvogel (1994)']);
  elseif (Vstretching == 2),
    disp(['Vstretching = ',num2str(Vstretching),                        ...
          '   Shchepetkin (2005)']);
  elseif (Vstretching == 3),
    disp(['Vstretching = ',num2str(Vstretching),                        ...
          '   Geyer (2009), BBL']);
  elseif (Vstretching == 4),
    disp(['Vstretching = ',num2str(Vstretching),                        ...
          '   Shchepetkin (2010)']);
  end
  if (kgrid == 1)
    disp(['   kgrid    = ',num2str(kgrid), '   at vertical W-points']);
  else
    disp(['   kgrid    = ',num2str(kgrid), '   at vertical RHO-points']);
  end
  disp(['   theta_s  = ',num2str(theta_s)]);
  disp(['   theta_b  = ',num2str(theta_b)]);
  disp(['   hc       = ',num2str(hc)]);
  disp(' ');
  disp(' S-coordinate curves: k, s(k), C(k)')
  disp(' ');


  if (kgrid == 1),
    for k=Nlev:-1:1,
      disp(['    ', ...
	    sprintf('%3g',k-1     ), '   ',                                 ...
	    sprintf('%20.12e',s(k)), '   ',                                 ...
	    sprintf('%20.12e',C(k))]);
    end,
  else
    for k=Nlev:-1:1,
      disp(['    ', ...
	    sprintf('%3g',k       ), '   ',                                 ...
	    sprintf('%20.12e',s(k)), '   ',                                 ...
	    sprintf('%20.12e',C(k))]);
    end,
  end,
  disp(' ');
end,

return

