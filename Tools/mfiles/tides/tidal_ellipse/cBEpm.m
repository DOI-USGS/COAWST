function [BEp, BEm]=cBEpm(g, f, sigma, nu, kappa, z, h)
%Evaluate the theoretical vertical profiles (or Bottom Ekman spiral profiles) 
%of tidal currents in the two rotary directions driven by half-unit of sea 
%surface gradients in the two directions respectively. Eddy viscosity is 
%assumed as vertically invariant. See tidal_ellipse.ps for more details.
%
%
%inputs:
%
%        g,  the gravity acceleration, 
%        f,  the Coriolis parameter,   
%       nu,  the eddy viscosity
%     kappa, the bottom frictional coefficient
%       z,   the vertical coordinates, can be a vector but must be
%            within [0 -h];
%       h,   the water depth, must be positive.
%
%       Note: except for z, all other inputs must be scalars.
%
%outputs:
%
%   BEp and BEm, the same dimensions of z,  the outputs for the vertical 
%                velocity profiles driven respectively by a unit of sea 
%                surface slope in the positive rotation direction and negative
%                rotation direction for when the eddy viscosity is vertically 
%                invariant. See the associated  document for more details.

if length(g)>1  | length(f)>1     | length(sigma)>1 | ...
   length(nu)>1 | length(kappa)>1 | length(h)>1
   error('inputs of g,f,sigma, nu, kappa, and h should be all scalars!')
end

if any(z/h>0) | any(z/h<-1)
   disp('z must be negative and must be within [0 -h]')
end

delta_e = sqrt(2*nu/f);  %Ekman depth
  alpha = (1+i)/delta_e*sqrt(1+sigma/f);
  beta  = (1+i)/delta_e*sqrt(1-sigma/f);

BEp = get_BE(g, alpha, h, z, nu, kappa);
BEm = get_BE(g, beta,  h, z,  nu, kappa);


%subfunction
%______________________________________________
function   BE=get_BE(g, alpha, h, z, nu, kappa)

       z = z(:);      
     z_h = z/h;     
      ah = alpha*h;
      az = alpha*z;
     ah2 = ah*2;
   anu_k = alpha*nu/kappa;    
   nu_kh = nu/(kappa*h);
     
    if abs(ah) < 1	%series solution
	 T = 10;
	 C = -g*h*h/(nu*(1+anu_k*tanh_v5_2(ah)))*2;   
	 A1 = (1-z_h.*z_h)/2+nu_kh;
	 B1 = exp(-ah)/(1+exp(-ah2)); 
	 B  = B1;
  	 series_sum=A1*B1;

	 for t = 2:T;      
		 t2=2*t;
		 A = (1-z_h.^t2)./t2+nu_kh;   
		 B = B*ah*ah/(t2-1)/(t2-2);
		 series_sum = series_sum+A*B;
	 end

	  BE = C*series_sum;

     else             %finite solution

          c = -g*h*h/nu;
          denom=(exp(az-ah)+exp(-(az+ah)))./(1+exp(-2*ah)); 
                                 % =cosh(az)/cosh(ah);
				 %but this a better way to evaluate it.

	  numer=1+anu_k*tanh_v5_2(ah);
	 % BE=c*(1-denom/numer);
          BE = c*((1-denom/numer)/(ah*ah));
     end
         
%end of subfunction
%
%Note tanh_v5_2 is a copy of tanh from Matlab v5.2, which has worked well!
%It seems that Matlab v5.3 has some bug(s) in tanh function! It cannot deal
%with large argument. try z=7.7249e02*(1+i), tanh(z) and tanh_v5_2(z) to 
%see the difference.

%Authorship Copyright:
%
%    The author of this program retains the copyright of this program, while
% you are welcome to use and distribute this program as long as you credit 
% the author properly and respect the program name itself. Particularly, 
% you are expected to retain the original author's name in this original 
% version of the program or any of its modified version that you might make.
% You are also expected not to essentially change the name of the programs 
% except for adding possible extension for your own version you might create, 
% e.g. ap2ep_xx is acceptable.  Any suggestions are welcome and enjoying my 
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

