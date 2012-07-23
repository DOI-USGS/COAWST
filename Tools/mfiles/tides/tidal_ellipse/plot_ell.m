function w=plot_ell(SEMA, ECC, INC, PHA, IND)
%
% An auxiliary function used in ap2ep and ep2ap for plotting
% tidal ellipse. The inputs, MA, ECC, INC and PHA are the output of 
% ap2ep and IND is a vector for indices for plotting a particular 
% ellipse, e.g., if IND=[2 3 1]; the ellipse corresponding to 
% the indices of (2,3,1) will be plotted.
%___________________________

 size_SEMA=size(SEMA);
  len_IND=length(IND);
  if IND   
       cmd=['sub2ind(size_SEMA' ];
       if len_IND==1   
	     titletxt=['Ellipse '];
       else
	     titletxt=['Ellipse ('];
       end

       for k=1:len_IND;
          cmd=[cmd ','num2str(IND(k))];
          if k<len_IND
         	titletxt=[titletxt num2str(IND(k)) ',']; 
	  elseif len_IND==1
                titletxt=[titletxt num2str(IND(k))]; 
          else
                titletxt=[titletxt num2str(IND(k)) ')']; 
          end
       end

       cmd=['n=' cmd ');'];
       eval(cmd);

       figure(gcf)
       clf
       do_the_plot(SEMA(n), ECC(n), INC(n), PHA(n));
       titletxt=[titletxt ',  (red) green (anti-) clockwise component'];
       title(titletxt);
   elseif len_IND
       msg1=['IND input contains zero element(s)!'];
       msg2=['No ellipse will be plotted.']; 
       disp(msg1);
       disp(msg2);
   end

%___________________________
%begin of plot subfunction
function w=do_the_plot(SEMA, ECC, INC, PHA)

   SEMI = SEMA.*ECC;
   Wp = (1+ECC)/2 .*SEMA;
   Wm = (1-ECC)/2 .*SEMA;
   THETAp = INC-PHA;
   THETAm = INC+PHA;

   %convert degrees into radians
   THETAp = THETAp/180*pi;
   THETAm = THETAm/180*pi;
   INC = INC/180*pi;
   PHA = PHA/180*pi;

   %Calculate wp and wm.
   wp = Wp.*exp(i*THETAp);
   wm = Wm.*exp(i*THETAm);
   
   dot = pi/18;
   ot = [0:dot:2*pi-dot];
   a = wp*exp(i*ot);
   b = wm*exp(-i*ot);   
   w = a+b;

   wmax = SEMA*exp(i*INC);
   wmin = SEMI*exp(i*(INC+pi/2));

   plot(real(w), imag(w))
   axis('equal');
   hold on
   plot([0 real(wmax)], [0 imag(wmax)], 'm');
   plot([0 real(wmin)], [0 imag(wmin)], 'm');
   xlabel('u');
   ylabel('v');
   plot(real(a), imag(a), 'r');
   plot(real(b), imag(b), 'g');
   hnd_a=line([0 real(a(1))], [0 imag(a(1))], 'color', 'r', 'marker','o');
   hnd_b=line([0 real(b(1))], [0 imag(b(1))], 'color', 'g', 'marker','o');
   hnd_w=line([0 real(w(1))], [0 imag(w(1))], 'color', 'b', 'marker','o');
   plot(real(a(1)), imag(a(1)), 'ro');
   plot(real(b(1)), imag(b(1)), 'go');
   plot(real(w(1)), imag(w(1)), 'bo');
   hnd_ab=line(real([a(1) a(1)+b(1)]), imag([a(1) a(1)+b(1)]), ... 
              'linestyle', '--', 'color', 'g');
   hnd_ba=line(real([b(1) a(1)+b(1)]), imag([b(1) a(1)+b(1)]), ... 
              'linestyle', '--', 'color', 'r');
  
   for n=1:length(ot);
         set(hnd_a, 'xdata', [0 real(a(n))], 'ydata', [0 imag(a(n))]);
         set(hnd_b, 'xdata', [0 real(b(n))], 'ydata', [0 imag(b(n))]);
         set(hnd_w, 'xdata', [0 real(w(n))], 'ydata', [0 imag(w(n))]);     
              hold on
         plot(real(a(n)), imag(a(n)), 'ro');
         plot(real(b(n)), imag(b(n)), 'go');
         plot(real(w(n)), imag(w(n)), 'bo');
         set(hnd_ab, 'xdata',real([a(n) a(n)+b(n)]), 'ydata', ... 
         imag([a(n) a(n)+b(n)]))
         set(hnd_ba, 'xdata',real([b(n) a(n)+b(n)]), 'ydata', ... 
         imag([b(n) a(n)+b(n)]))
   end

   hold off

%end of plot subfunction
%---------------------------
%
%

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

