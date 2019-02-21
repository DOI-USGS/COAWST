function [Sf,Sfth,Freq,Dir] = WWIIIpartitioned2full2D(Hsig,Tp,theta,sigma,gamma,freq,dth);
%
% function [Sf,Sfth,Freq,Dir] = WWIIIpartitioned2full2D(Hsig,Tp,theta,sigma,gamma,freq);
% Usage: 
% 
% This function takes bulk wave parameters (Significant Wave Ht., Peak Period, 
% Peak Direction and Directional Spreading) to create a frequency-directional ocean wave spectrum
% using a combination of Jonswap spectrum and a cos^2s(theta) distribution.
%
% Developed because of the need to convert the partitioned wave information  by WW III to a full 2-D spectra
% According to this method a JONSWAP type spectrum is created for each partition and then a directional 
% distribution is applied to it. Finally all spectra are added together ensuring energy preservation.
%
%  INPUT: 
% 		Hsig   = Significant wave height (in m)
% 		Tp     = Peak Period (in sec)
% 		theta  = Peak Direction (in Degrees)
% 		sigma  = Spectral width/Directional spreading in degrees
% 		gammma = [Optional, default value =3.3] Peakedness Parameter (Decides the peakedness of a Jonswap spectrum, values from 1-7)
%       freq   = [Optional] Frequency array for the full spectrum. If no frequency array is provided
%               the function uses the default one from 0.01 to 1Hz with a df=0.01Hz
%       dth    = [Optional] Directional width for the full spectrum (in degrees). If no
%                value is provided, the function uses the default value of 6 degrees
%  OUTPUT:
% 		Sf     = Spectral density as a function of frequency (m^2/Hz)
% 		Sfth   = Frequency Directional Wave Spectrum (m^2/Hz/Deg)
% 		Freq   = Frequency Range (in Hz)
% 		Dir    = Direction Range (in Degrees), with radial resolution dth
% 		         in degrees
% 
%  EXAMPLES:
%
%      [Sf,Sfth,Freq,Dir] = WWIIIpartitioned2full2D(1,5,0,30,[],[],[]);
%      [Sf,Sfth,Freq,Dir] = WWIIIpartitioned2full2D(1,5,0,30,5,[],[]);
%      [Sf,Sfth,Freq,Dir] = WWIIIpartitioned2full2D([1,2,3],[5 7 10],[0 120,30],[30 22 40]);
%      [Sf,Sfth,Freq,Dir] = WWIIIpartitioned2full2D(1,5,0,30,[],[0.01:0.01:0.5]);
%  USES Internal functions (i) create_jonswap_spectrum.m, (ii) Directional_Distribution.m
%
% Nirnimesh Kumar and George Voulgaris
% Coastal Processes and Sediment Dynamics Lab
% Dept. of Earth and Ocean Sciences,
% Univ. of South Carolina, Columbia, SC
% 03/06/2013

% NKumar, 10/21/2014
% (a) Cleaned code and changed directional distribution following Kuik et al. (1988)
% (b) Directional spread power (i.e., s) needs to be an integer to get 
% positive values of directional distribution function, which requires
% rounding of the value of s, as done in line number 167

if (nargin<7 || isempty(dth)==1)
    dth = 6; %(This is same as the dtheta I would use for SWAN)  
end

if (nargin<6 || isempty(freq)==1)
    freq = 0.01:0.01:1;
end

if (nargin<5 || isempty(gamma)==1)
    gamma= 3.3;
end

if nargin<4
	error('Function needs Significant Wave Ht, Peak Period, Peak Direction,Directional Spreading and Gamma Parameter')
end

N = length(Hsig);
%
for i=1:N
    [Sfo,Freq]   = create_jonswap_spectrum(Hsig(i),Tp(i),gamma,freq);
    [Sftho,Dir]  = Directional_Distribution(Sfo,theta(i),sigma(i),dth);
    if i==1;
        Sfth =zeros(size(Sftho));
        Sf   =zeros(size(Sfo));
    end
    Sfth         = Sfth+Sftho;
    Sf           = Sf+Sfo;
end
end

function [Sf,freq]=create_jonswap_spectrum(Hsig,Tp,gamma,freq);
%
%  function [Sf]=create_jonswap_spectrum(Hsig,Tp,gamma);
%  Usage: This function creates a Jonswap Spectrum on the basis of given
%  significant wave height (Hsig), peak wave period (Tp) and peakedness parameter (gamma).
%  Another possibility is to use the WAFO Toolbox (http://www.maths.lth.se/matstat/wafo/) which 
%  provides a more accurate implementation of JONSWAP spectrum:
%
%  INPUT:
%		Hsig = Significant wave height (in m)
%		Tp   = Peak wave period (in sec)
%		gamma= Peakedness parameter (default value 3.3).
%		freq = Frequency distribution (default 0:0.01:1)
%
%	OUTPUT:
%		Sf   = Spectral density as a function of frequency (m^2/Hz)
%
% Nirnimesh Kumar and George Voulgaris
% Coastal Processes and Sediment Dynamics Lab
% Dept. of Earth and Ocean Sciences,
% Univ. of South Carolina, Columbia, SC
% 03/06/2013
%
if nargin<4
    freq = 0.01:0.01:1;
end
if nargin<3
	error('Function needs at least Significant Wave Ht, Peak Period and Gamma Parameter')
end

g     = 9.81;           % Acceleration due to gravity
Hsig  = Hsig(:);        % Significant wave height
Tp    = Tp(:);          % Peak Period
gamma = gamma(:);       % Peakedness parameter, gamma
freq  = freq(:);        
omegap= 2*pi./Tp;       % Peak angular frequency
omega = 2*pi*freq;      % Angular frequency
domega= diff(omega);
domega= [domega(1);domega];

Sw     = zeros(length(freq),1); %Create array for spectral density
Sf     = zeros(length(freq),1);
sigma  = zeros(length(freq),1);

j       = omega>omegap;
sigma(j)= 0.09;
j       = omega<=omegap;
sigma(j)= 0.07;

a       = exp(-((omega-omegap).^2)./(2*(omegap.^2).*(sigma.^2)));
beta    = 5/4;
Sw      = (1./(omega.^5)).*(exp(-beta*(omegap.^4)./(omega.^4))).*(gamma.^a);
NormFac = ((Hsig./g).^2)./16./sum(Sw.*domega);

%Normalized spectrum:
Sw      = Sw*NormFac*g^2; 
Sf      = 2*pi*Sw;           % Units are m^2/Hz

end

function [Sfth,Dir]=Directional_Distribution(Sf,theta,sigma,dth) 
%
% function [Sfth]=Directional_Distribution(Sf,theta,sigma) 
%
% 	Usage: 
% This function considers the frequency spectrum, peak direction and
% directional spreading (defined as spectral width in Wavewatch Data) and
% applies a cos^2s(theta-theta_peak) distribution where theta varies from 0
% to 360 degrees.
%
%	INPUT:
%		Sf    = Spectral Density as a function of frequency (units m^2/Hz)
%		theta = Peak Direction (of a partition) in degrees
%		sigma = Spectral Width/Directional Spreading in degrees
%       dth   = radial resolution of direction in degs
%	OUTPUT:
%		Sfth  = Frequency Directional Wave Spectrum (m^2/Hz/Deg)
%
% Nirnimesh Kumar and George Voulgaris
% Coastal Processes and Sediment Dynamics Lab
% Dept. of Earth and Ocean Sciences,
% Univ. of South Carolina, Columbia, SC
% 03/06/2013
%

Sf  = Sf(:);
Th  = (pi/180)*[dth:dth:360]';           % direction from 0 to 360 degrees
s   = round((2./((pi*sigma/180).^2))-1); %See Kuik et al. (1988) 
if s>80
    s=80;
end
if s<1
    s=1;
end
aa   = gamma(s+1);
bb   = gamma(2*s+1);
Fth  = (2^(2*s-1)/pi)*(aa^2/bb)*cos((Th-theta*pi/180)/2).^(2*s);
Fth  = real(Fth*pi/180); %This makes the units consistent with m^2/Hz/Deg (For e.g., sum(Fth*dth) must be 1)
Sfth = (Sf*Fth');
Dir  = Th*180/pi;
end

%%%%%%%%%%%%%%%%%%%%%%Old directional distributions%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convert directional spreading to spreading factor expnent (s) empirically
% (This equation is based on relation shown in SWAN Manual, also see Kuik et al., 1988)
%
% Dum  = -1.2578*(log10(sigma)).^5 + 4.8416*(log10(sigma)).^4 - 7.0731*(log10(sigma)).^3 ...
%       +4.7358*(log10(sigma)).^2 - 3.4120*log10(sigma) + 3.6599;
% 
% s    = max(1,round(0.5*(10.^(Dum))));
%s    = 5*max(1,round(0.5*(10.^(Dum))));
% Dum   = - 0.89333*(log10(sigma)).^4 + 4.2303*(log10(sigma)).^3 - 7.5468*(log10(sigma)).^2 + 3.9694*(log10(sigma)) + 2.0547;
% s     = 10.^Dum;