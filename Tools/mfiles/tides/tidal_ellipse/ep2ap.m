function [Au, PHIu, Av, PHIv, w]=ep2ap(SEMA, ECC, INC, PHA, plot_demo)
% Convert tidal ellipse parameters into amplitude and phase lag parameters.
% Its inverse is app2ep.m. Please refer to app2ep for the meaning of the 
% inputs and outputs.
%
% Zhigang Xu
% Oct. 20, 2000
%
% Document:  tidal_ellipse.ps
% 
if nargin < 5
     plot_demo=0;  % by default, no plot for the ellipse
end

   Wp = (1+ECC)/2 .*SEMA;
   Wm = (1-ECC)/2 .*SEMA;
   THETAp = INC-PHA;
   THETAm = INC+PHA;

   %convert degrees into radians
   THETAp = THETAp/180*pi;
   THETAm = THETAm/180*pi;

   %Calculate wp and wm.
   wp = Wp.*exp(i*THETAp);
   wm = Wm.*exp(i*THETAm);
   
   if nargout == 5
      ndot=36;
      dot = 2*pi/ndot;
      ot = [0:dot:2*pi-dot];
      w = wp(:)*exp(i*ot)+wm(:)*exp(-i*ot);
      w=reshape(w, [size(wp) ndot]);
   end

   % Calculate cAu, cAv --- complex amplitude of u and v
   cAu = wp+conj(wm);
   cAv = -i*(wp-conj(wm));
   Au  = abs(cAu);
   Av  = abs(cAv);   
   PHIu = -angle(cAu)*180/pi;
   PHIv = -angle(cAv)*180/pi;
   
   % flip angles in the range of [-180 0) to the range of [180 360).
   id = PHIu < 0; PHIu(id) = PHIu(id) + 360;
   id = PHIv < 0; PHIv(id) = PHIv(id) + 360;

  if any(plot_demo)   
     plot_ell(SEMA,ECC,INC,PHA,plot_demo)
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


