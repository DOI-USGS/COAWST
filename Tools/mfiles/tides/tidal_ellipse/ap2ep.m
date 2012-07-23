function [SEMA,  ECC, INC, PHA, w]=ap2ep(Au, PHIu, Av, PHIv, plot_demo)
% Convert tidal amplitude and phase lag (ap-) parameters into tidal ellipse
% (e-) parameters. Please refer to ep2app for its inverse function.
% 
% Usage:
%
% [SEMA,  ECC, INC, PHA, w]=app2ep(Au, PHIu, Av, PHIv, plot_demo)
%
% where:
%
%     Au, PHIu, Av, PHIv are the amplitudes and phase lags (in degrees) of 
%     u- and v- tidal current components. They can be vectors or 
%     matrices or multidimensional arrays.
%            
%     plot_demo is an optional argument, when it is supplied as an array 
%     of indices, say [i j k l], the program will plot an  ellipse 
%     corresponding to Au(i,j, k, l), PHIu(i,j,k,l), Av(i,j,k,l), and 
%     PHIv(i,j,k,l); 
%     
%     Any number of dimensions are allowed as long as your computer 
%     resource can handle.     
%     
%     SEMA: Semi-major axes, or the maximum speed;
%     ECC:  Eccentricity, the ratio of semi-minor axis over 
%           the semi-major axis; its negative value indicates that the ellipse
%           is traversed in clockwise direction.           
%     INC:  Inclination, the angles (in degrees) between the semi-major 
%           axes and u-axis.                        
%     PHA:  Phase angles, the time (in angles and in degrees) when the 
%           tidal currents reach their maximum speeds,  (i.e. 
%           PHA=omega*tmax).
%          
%           These four e-parameters will have the same dimensionality 
%           (i.e., vectors, or matrices) as the input ap-parameters. 
%
%     w:    Optional. If it is requested, it will be output as matrices
%           whose rows allow for plotting ellipses and whose columns are  
%           for different ellipses corresponding columnwise to SEMA. For
%           example, plot(real(w(1,:)), imag(w(1,:))) will let you see 
%           the first ellipse. You may need to use squeeze function when
%           w is a more than two dimensional array. See example.m. 
%
% Document:   tidal_ellipse.ps

if nargin < 5
     plot_demo=0;  % by default, no plot for the ellipse
end


% Assume the input phase lags are in degrees and convert them in radians.
   PHIu = PHIu/180*pi;
   PHIv = PHIv/180*pi;

% Make complex amplitudes for u and v
   i = sqrt(-1);
   u = Au.*exp(-i*PHIu);
   v = Av.*exp(-i*PHIv);

% Calculate complex radius of anticlockwise and clockwise circles:
   wp = (u+i*v)/2;      % for anticlockwise circles
   wm = conj(u-i*v)/2;  % for clockwise circles
% and their amplitudes and angles
   Wp = abs(wp);
   Wm = abs(wm);
   THETAp = angle(wp);
   THETAm = angle(wm);
   
% calculate e-parameters (ellipse parameters)
    SEMA = Wp+Wm;              % Semi  Major Axis, or maximum speed
    SEMI = Wp-Wm;              % Semin Minor Axis, or minimum speed
     ECC = SEMI./SEMA;          % Eccentricity

    PHA = (THETAm-THETAp)/2;   % Phase angle, the time (in angle) when 
                               % the velocity reaches the maximum
    INC = (THETAm+THETAp)/2;   % Inclination, the angle between the 
                               % semi major axis and x-axis (or u-axis).

    % convert to degrees for output
    PHA = PHA/pi*180;         
    INC = INC/pi*180;         
    THETAp = THETAp/pi*180;
    THETAm = THETAm/pi*180;
    
  % flip THETAp and THETAm, PHA, and INC in the range of 
  % [-pi, 0) to [pi, 2*pi), which at least is my convention.
    id = THETAp < 0;   THETAp(id) = THETAp(id)+360;
    id = THETAm < 0;   THETAm(id) = THETAm(id)+360;
    id = PHA < 0;      PHA(id) = PHA(id)+360;
    id = INC < 0;      INC(id) = INC(id)+360;
    

  if nargout == 5
     ndot=36;
     dot=2*pi/ndot;
     ot=[0:dot:2*pi-dot];
     w=wp(:)*exp(i*ot)+wm(:)*exp(-i*ot);
     w=reshape(w, [size(Au) ndot]);
  end


 if any(plot_demo)
    plot_ell(SEMA, ECC, INC, PHA, plot_demo)
 end

%Authorship Copyright:
%
%    The author of this program retains the copyright of this program, while
% you are welcome to use and distribute this program as long as you credit 
% the author properly and respect the program name itself. Particularly, 
% you are expected to retain the original author's name in this original 
% version of the program or any of its modified version that you might make.
% You are also expected not to essentially change the name of the programs 
% except for adding possible extension for your own version you might create, 
% e.g. app2ep_xx is acceptable.  Any suggestions are welcome and enjoying my 
% program(s)!
%
%
%Author Info:
%_______________________________________________________________________
%  Zhigang Xu, Ph.D.                            
%  (pronounced as Tsi Gahng Hsu)
%  Research Scientist
%  Coastal Circulation                   
%  Bedford Institute of Oceanography     
%  1 Challenge Dr.
%  P.O. Box 1006                    Phone  (902) 426-2307 (o)       
%  Dartmouth, Nova Scotia           Fax    (902) 426-7827            
%  CANADA B2Y 4A2                   email zhigangx@emerald.bio.dfo.ca    
%                                         zhigang_xu_98@yahoo.com
%_______________________________________________________________________
%
%Release Date: Nov. 2000


