function [astro,ader] = t_astron(jd)
% T_ASTRON Computes astronomical Variables
% [A,ADER] = ASTRON(JD) computes the astronomical variables 
%            A=[tau,s,h,p,np,pp] (cycles) 
%  and their time derivatives 
%            ADER=[dtau,ds,dh,dp,dnp,dpp] (cycles/day) 
%  at the matlab time JD (UTC, but see code for details) where
%
%	tau = lunar time
%	s = mean longitude of the moon
%	h = mean longitude of the sun
%	p = mean longitude of the lunar perigee 
%	np = negative of the longitude of the mean ascending node
%	pp = mean longitude of the perihelion (solar perigee)   
%

%
%    The formulae for calculating these ephemerides (other than tau) 
%    were taken from pages 98 and 107 of the Explanatory Supplement to
%    the Astronomical Ephemeris and the American Ephemeris and Nautical 
%    Almanac (1961). They require EPHEMERIS TIME (ET), now TERRESTRIAL 
%    TIME (TT) and are based on observations made in the 1700/1800s.
%    In a bizarre twist, the current definition of time is derived
%    by reducing observations of planetary motions using these formulas.
%
%    The current world master clock is INTERNATIONAL ATOMIC TIME (TAI).
%    The length of the second is based on inverting the actual 
%    locations of the planets over the period 1956-65 into "time" 
%    using these formulas, and an offset added to keep the scale 
%    continuous with previous defns. Thus
%
%                     TT = TAI + 32.184 seconds.
%
%    Universal Time UT is a time scale that is 00:00 at midnight (i.e.,
%    based on the earth's rotation rather than on planetary motions).
%    Coordinated Universal Time (UTC) is kept by atomic clocks, the 
%    length of the second is the same as for TAI but leap seconds are
%    inserted at intervals so that it provides UT to within 1 second. 
%    This is necessary because the period of the earth's rotation is 
%    slowly increasing (the day was exactly 86400 seconds around 1820, 
%    it is now about 2 ms longer). 22 leap seconds have been added in 
%    the last 27 years.
%
%    As of 1/1/99,    TAI = UTC + 32 seconds.
%       
%    Thus,             TT = UTC + 62.184 seconds
%
%    GPS time was synchronized with UTC 6/1/1980 ( = TAI - 19 secs), 
%    but is NOT adjusted for leap seconds. Your receiver might do this
%    automatically...or it might not.
%
%    Does any of this matter? The moon longitude is the fastest changing
%    parameter at 13 deg/day. A time error of one minute implies a
%    position error of less than 0.01 deg. This would almost always be 
%    unimportant for tidal work.
%
%    The lunar time (tau) calculation requires UT as a base.  UTC is 
%    close enough - an error of 1 second, the biggest difference that
%    can occur between UT and UTC, implies a Greenwich phase error of 
%    0.01 deg.  In Doodson's definition (Proc R. Soc. A, vol 100, 
%    reprinted in International Hydrographic Review, Appendix to 
%    Circular Letter 4-H, 1954) mean lunar time is taken to begin at 
%    "lunar midnight". 

% B. Beardsley  12/29/98, 1/11/98
% R. Pawlowicz  9/1/01
% Version 1.0


% Compute number of days from epoch of 12:00 UT Dec 31, 1899.
% (January 0.5 1900 ET)
d=jd(:)'-datenum(1899,12,31,12,0,0);
D=d/10000;

% Compute astronomical constants at time d1.
args=[ones(size(jd));
      d;
      D.*D;
      D.^3];

% These are the coefficients of the formulas in the Explan. Suppl.

sc= [ 270.434164,13.1763965268,-0.0000850, 0.000000039];
hc= [ 279.696678, 0.9856473354, 0.00002267,0.000000000];
pc= [ 334.329556, 0.1114040803,-0.0007739,-0.00000026];
npc=[-259.183275, 0.0529539222,-0.0001557,-0.000000050];
%  first coeff was 281.220833 in Foreman but Expl. Suppl. has 44.
ppc=[ 281.220844, 0.0000470684, 0.0000339, 0.000000070];

% Compute the parameters; we only need the factional part of the cycle.
astro=rem( [sc;hc;pc;npc;ppc]*args./360.0 ,1);

% Compute lunar time tau, based on fractional part of solar day.
% We add the hour angle to the longitude of the sun and subtract the
% longitude of the moon.
tau=rem(jd(:)',1)+astro(2,:)-astro(1,:);
astro=[tau;astro];

% Compute rates of change.
dargs=[zeros(size(jd));
       ones(size(jd));
       2.0e-4.*D;
       3.0e-4.*D.*D];

ader=[sc;hc;pc;npc;ppc]*dargs./360.0;

dtau=1.0+ader(2,:)-ader(1,:);

ader=[dtau;ader];

