function [name,freq,amp]=t_equilib(lat);
% T_EQUILIB Equilibrium amplitude of the tidal potential
% [NAME,FREQ,AMPLITUDE]=T_EQUILIB(LAT) returns vectors with the 
% NAME of tidal constituents, their FREQ (in cph), and their 
% equilibrium AMPLITUDE in the tidal potential as a function of 
% LATitude (degrees). If LAT is a vector, then AMPLITUDE is a 
% matrix in which each column corresponds to a specific latitude.
%
% If no output arguments are specified, the equilibrium spectrum
% is plotted.

% R. Pawlowicz 9/11/99
% Version 1.0


const=t_get18consts;

g=9.81;            % m/s^2;
erad=6365;         % km
earthmoond=3.84e5; % km
Mmoon=7.38e22;     % kg
Mearth=5.977e24;   % kg
Gravconst=6.658e-11;  % m^3/kg/s^2

% There appears to be a typo in Godin's text, and this
% should likely be checked against Doodson's original.
% This is what I *think* it should be.
G=3/4*Mmoon*(erad*1e3)^3/(earthmoond*1e3)^2/Mearth;

% The 4/3 is to correct for the 3/4 in G
gfac=Gravconst*Mearth/(erad*1e3)^2*(4/3);

jk=finite(const.doodsonamp);


freq=const.freq(jk);
name=const.name(jk,:);
	
slat=sin(lat(:)'*pi/180);
clat=cos(lat(:)'*pi/180);

G1=zeros(6,length(clat));

% Latitude dependence of amplitude for various species -
% + for A, -for B (from Godin, 1972).

G1(3+0,:)=    0.5*G*(1-3*slat.^2);
G1(3-1,:)=      2*G*slat.*clat;
G1(3+1,:)= .72618*G*clat.*(1-5*slat.^2);
G1(3-2,:)=2.59808*G*slat.*clat.^2;
G1(3+2,:)=        G*clat.^2;
G1(3+3,:)=        G*clat.^3;

	
amp=abs(const.doodsonamp(jk,ones(1,length(clat)))/gfac.*G1(const.doodsonspecies(jk)+3,:));
 
 
if nargout==0,

 cols=[1 0 0;
       0 1 0;
       0 0 1;
       .5 0 0;
       0 .5 0;
       0 0 .5];
       
 set(gcf,'defaultaxescolororder',cols(const.doodsonspecies(jk)+3,:));
 semilogy(24*[freq,freq]',[repmat(min(amp),length(amp),1) amp]');
 
 cnam=cellstr(name);
 for k=1:length(cnam),
   cnam{k}=deblank(cnam{k});
   ff=min([find(abs(cnam{k}(2:end))>=abs('0') & abs(cnam{k}(2:end))<=abs('9'))+1,length(cnam{k})+1]);
   cnam{k}=[ cnam{k}(1:ff-1) '_{' cnam{k}(ff:end) '}'];
 end;
 text(freq*24,amp,cnam,'vertical','bottom','horiz','center','fontangle','italic','fontweight','bold',...
       'clip','on','fontsize',9);
 xlabel('Frequency (cpd)');
 ylabel('Potential (m)');
 set(gca,'tickdir','out','ylim',[min(amp) max(amp)*2]);
 
end;






    
