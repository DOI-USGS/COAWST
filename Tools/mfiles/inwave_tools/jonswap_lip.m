function [ S, Amp, Phase ] = jonswap( Ohm, Hs, Tp)
%JONSWAP - Calculates the wave spectrum values for a JONSWAP spectrum
%
% For the Lip 2e case I used:
% f   = [0.0001:0.01:1.0];
% Ohm = 2*pi./f;
% Hs  = 1.4;
% Tp  = 5;

 jcw=1;

 wp = 2*pi/Tp;
 Gamma = 3.3;
 for x = 1:length(Ohm)
     if Ohm(x)<wp
         Sigma = 0.07;
     else
         Sigma = 0.09;
     end
     ff=2*pi/Ohm(x);
     fp=1/Tp;
     A = exp(-((Ohm(x)/wp-1)/(Sigma*sqrt(2)))^2);
     if (jcw)
       S(x) = 0.0081*9.81^2/((2*3.141)^4)/(ff^5)*exp(-5/4*(fp^4)/(ff^4))*Gamma^A;
     else
       S(x) = 320*Hs^2*Ohm(x)^-5/Tp^4*exp(-1950*Ohm(x)^-4/Tp^4)*Gamma^A;
     end
 end

 % Determine the frequency step from the frequency vector. Note that the
 % highest frequency step is extrapolated.
 domg = zeros( size(Ohm) );
 domg(1:end-1) = diff( Ohm );
 domg(end) = domg(end-1);

 % Determine the amplitudes from the spectral values
 Amp = sqrt( 2 * S .* domg );

 % Random phases
 Phase = rand(1,length(Ohm))*2*pi;

 if (jcw)
   % scale S so the integral of S*df gives you Hmo
   f=2*pi./Ohm;
   m0=trapz(f,S);
   Hact1=4.004*sqrt(m0)
   Scale=Hs^2/16/m0;
   S=S.*Scale;
   m0=trapz(f,S);
   Hend=4.004*sqrt(m0)
 end

end





