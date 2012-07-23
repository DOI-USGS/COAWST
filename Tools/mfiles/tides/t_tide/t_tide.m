function [nameu,fu,tidecon,xout]=t_tide(xin,varargin);
% T_TIDE Harmonic analysis of a time series
% [NAME,FREQ,TIDECON,XOUT]=T_TIDE(XIN) computes the tidal analysis 
% of the (possibly complex) time series XIN.
%
% [TIDESTRUC,XOUT]=T_TIDE(XIN) returns the analysis information in
% a structure formed of NAME, FREQ, and TIDECON.
%
% XIN can be scalar (e.g. for elevations), or complex ( =U+sqrt(-1)*V
% for eastward velocity U and northward velocity V.
%
% Further inputs are optional, and are specified as property/value pairs
% [...]=T_TIDE(XIN,property,value,property,value,...,etc.)
%      
% These properties are:
%
%       'interval'       Sampling interval (hours), default = 1. 
%          
%   The next two are required if nodal corrections are to be computed,
%   otherwise not necessary. If they are not included then the reported
%   phases are raw constituent phases at the central time. 
%
%   If your time series is longer than 18.6 years then nodal corrections
%   are not made -instead we fit directly to all satellites (start time
%   is then just used to generate Greenwich phases).
%
%       'start time'     [year,month,day,hour,min,sec]
%                        - min,sec are optional OR 
%                        decimal day (matlab DATENUM scalar)
%       'latitude'       decimal degrees (+north) (default: none).
%
%   Where to send the output.
%       'output'         where to send printed output:
%                        'none'    (no printed output)
%                        'screen'  (to screen) - default
%                        FILENAME   (to a file)
%
%   Correction factor for prefiltering.
%       'prefilt'        FS,CORR
%                        If the time series has been passed through
%                        a pre-filter of some kind (say, to reduce the
%                        low-frequency variability), then the analyzed
%                        constituents will have to be corrected for 
%                        this. The correction transfer function 
%                        (1/filter transfer function) has (possibly 
%                        complex) magnitude CORR at frequency FS (cph). 
%                        Corrections of more than a factor of 100 are 
%                        not applied; it is assumed these refer to tidal
%                        constituents that were intentionally filtered 
%                        out, e.g., the fortnightly components.
%
%   Adjustment for long-term behavior ("secular" behavior).
%       'secular'        'mean'   - assume constant offset (default).
%                        'linear' - get linear trend.
%                     
%   Inference of constituents.
%       'inference'      NAME,REFERENCE,AMPRAT,PHASE_OFFSET
%                        where NAME is an array of the names of 
%                        constituents to be inferred, REFERENCE is an 
%                        array of the names of references, and AMPRAT 
%                        and PHASE_OFFSET are the amplitude factor and
%                        phase offset (in degrees)from the references. 
%                        NAME and REFERENCE are Nx4 (max 4 characters
%                        in name), and AMPRAT and PHASE_OFFSET are Nx1
%                        (for scalar time series) and Nx2 for vector 
%                        time series (column 1 is for + frequencies and
%                        column 2 for - frequencies).
%                        NB - you can only infer ONE unknown constituent
%                        per known constituent (i.e. REFERENCE must not 
%                        contain multiple instances of the same name).
%
%   Shallow water constituents
%       'shallow'        NAME
%                        A matrix whose rows contain the names of 
%                        shallow-water constituents to analyze.
%
%   Resolution criterions for least-squares fit.        
%       'rayleigh'       scalar - Rayleigh criteria, default = 1.
%                        Matrix of strings - names of constituents to
%                                   use (useful for testing purposes).
%  
%   Calculation of confidence limits.
%       'error'          'wboot'  - Boostrapped confidence intervals 
%                                   based on a correlated bivariate 
%                                   white-noise model.
%                        'cboot'  - Boostrapped confidence intervals 
%                                   based on an uncorrelated bivariate 
%                                   coloured-noise model (default).
%                        'linear' - Linearized error analysis that 
%                                   assumes an uncorrelated bivariate 
%                                   coloured noise model. 
%                                   
%   Computation of "predicted" tide (passed to t_predic, but note that
%                                    the default value is different).
%       'synthesis'      0 - use all selected constituents
%                        scalar>0 - use only those constituents with a 
%                                   SNR greater than that given (1 or 2 
%                                   are good choices, 2 is the default).
%                              <0 - return result of least-squares fit 
%                                   (should be the same as using '0', 
%                                   except that NaN-holes in original 
%                                   time series will remain and mean/trend
%                                   are included).
%
%   Least squares soln computational efficiency parameter
%	'lsq'		'direct'  - use A\x fit
%			'normal'  - use (A'A)\(A'x) (may be necessary
%				    for very large input vectors since
%                                   A'A is much smaller than A)
%			'best'	  - automatically choose based on
%				    length of series (default).
%
%       It is possible to call t_tide without using property names,
%       in which case the assumed calling sequence is
%
%          T_TIDE(XIN,INTERVAL,START_TIME,LATITUDE,RAYLEIGH)
%
%
%  OUTPUT: 
%
%    nameu=list of constituents used
%    fu=frequency of tidal constituents (cycles/hr)
%    tidecon=[fmaj,emaj,fmin,emin,finc,einc,pha,epha] for vector xin
%           =[fmaj,emaj,pha,epha] for scalar (real) xin
%       fmaj,fmin - constituent major and minor axes (same units as xin)       
%       emaj,emin - 95% confidence intervals for fmaj,fmin
%       finc - ellipse orientations (degrees)
%       einc - 95% confidence intervals for finc
%       pha - constituent phases (degrees relative to Greenwich)
%       epha - 95% confidence intervals for pha
%    xout=tidal prediction
%
% Note: Although missing data can be handled with NaN, it is wise not
%       to have too many of them. If your time series has a lot of 
%       missing data at the beginning and/or end, then truncate the 
%       input time series.  The Rayleigh criterion is applied to 
%       frequency intervals calculated as the inverse of the input 
%       series length.
%
% A description of the theoretical basis of the analysis and some
% implementation details can be found in:
%
% Pawlowicz, R., B. Beardsley, and S. Lentz, "Classical Tidal 
%   "Harmonic Analysis Including Error Estimates in MATLAB 
%    using T_TIDE", Computers and Geosciences, 28, 929-937 (2002).
%
% (citation of this article would be appreciated if you find the
%  toolbox useful).


% R. Pawlowicz 11/8/99 - Completely rewritten from the transliterated-
%                        to-matlab IOS/Foreman fortran code by S. Lentz
%                        and B. Beardsley.
%              3/3/00  - Redid errors to take into account covariances 
%                        between u and v errors.
%              7/21/00 - Found that annoying bug in error calc! 
%              11/1/00 - Added linear error analysis.
%              8/29/01 - Made synth=1 default, also changed behavior 
%                        when no lat/time given so that phases are raw
%                        at central time. 
%              9/1/01  - Moved some SNR code to t_predic.
%              9/28/01 - made sure you can't choose Z0 as constituent.
%              6/12/01 - better explanation for variance calcs, fixed
%                        bug in typed output (thanks Mike Cook).
%              8/2/03 - Added block processing for long time series (thanks
%                       to Derek Goring).
%              9/2/03 - Beta version of 18.6 year series handling
%              12/2/03 - Bug (x should be xin) fixed thanks to Mike Cook (again!)
%              4/3/11 - Changed (old) psd to (new) pwelch calls, also
%                       isfinite for finite.
%              23/3/11 - Corrected my conversion from psd to pwelch, thanks
%                       to Dan Codiga and (especially) Evan Haug!

%
% Version 1.3



% ----------------------Parse inputs-----------------------------------

ray=1;
dt=1;
fid=1;
stime=[];
lat=[];
corr_fs=[0 1e6];
corr_fac=[1  1];
secular='mean';
inf.iname=[];
inf.irefname=[];
shallownames=[];
constitnames=[];
errcalc='cboot';
synth=2;
lsq='best';

k=1;
while length(varargin)>0,
  if ischar(varargin{1}),
    switch lower(varargin{1}(1:3)),
      case 'int',
        dt=varargin{2};
      case 'sta',
        stime=varargin{2};
	if length(stime)>1, 
	  stime=[stime(:)' zeros(1,6-length(stime))]; 
	  stime=datenum(stime(1),stime(2),stime(3),stime(4),stime(5),stime(6));
	end;
      case 'lat',
         lat=varargin{2};
      case 'out',
         filen=varargin{2};
	 switch filen,
	   case 'none',
	     fid=-1;
	   case 'screen',
	     fid=1;
	   otherwise
	     [fid,mesg]=fopen(filen,'w');
	     if fid==-1, error(msg); end;
	  end;
      case 'ray',
         if isnumeric(varargin{2}),
           ray=varargin{2};
	 else
	   constitnames=varargin{2};
	   if iscellstr(constitnames), constitnames=char(constitnames); end;
	 end;
       case 'pre',
         corr_fs=varargin{2};
	 corr_fac=varargin{3};
         varargin(1)=[];
      case 'sec',
         secular=varargin{2};
      case 'inf',
         inf.iname=varargin{2};
	 inf.irefname=varargin{3};
	 inf.amprat=varargin{4};
	 inf.ph=varargin{5};
	 varargin(1:3)=[];
      case 'sha',
         shallownames=varargin{2};
      case 'err',
         errcalc=varargin{2};
      case 'syn',
         synth=varargin{2};
      case 'lsq',
         lsq=varargin{2};	 
      otherwise,
         error(['Can''t understand property:' varargin{1}]);
    end;
    varargin([1 2])=[]; 
  else  
    switch k,
      case 1,
        dt=varargin{1};
      case 2,
        stime=varargin{1};
      case 3,
        lat=varargin{1};
      case 4,
        ray=varargin{1};
      otherwise
        error('Too many input parameters');
     end;
     varargin(1)=[];
  end;
  k=k+1;
end;
 
[inn,inm]=size(xin);
if ~(inn==1 | inm==1), error('Input time series is not a vector'); end;

xin=xin(:); % makes xin a column vector
nobs=length(xin);

if strcmp(lsq(1:3),'bes'),  % Set matrix method if auto-choice.
 if nobs>10000,
    lsq='normal';
 else
    lsq='direct';
 end;
end;
 
if nobs*dt> 18.6*365.25*24,  % Long time series
  longseries=1; ltype='full';
else
  longseries=0; ltype='nodal';
end;
        		
nobsu=nobs-rem(nobs-1,2);% makes series odd to give a center point

t=dt*([1:nobs]'-ceil(nobsu/2));  % Time vector for entire time series,
                                 % centered at series midpoint. 

if ~isempty(stime),
  centraltime=stime+floor(nobsu./2)./24.0*dt;
else
  centraltime=[];
end;

% -------Get the frequencies to use in the harmonic analysis-----------

[nameu,fu,ju,namei,fi,jinf,jref]=constituents(ray/(dt*nobsu),constitnames,...
                                           shallownames,inf.iname,inf.irefname,centraltime);

mu=length(fu); % # base frequencies
mi=length(fi); % # inferred

% Find the good data points (here I assume that in a complex time 
% series, if u is bad, so is v).

gd=find(isfinite(xin(1:nobsu)));
ngood=length(gd);
fprintf('   Points used: %d of %d\n',ngood,nobs)



%----------------------------------------------------------------------
% Now solve for the secular trend plus the analysis. Instead of solving
% for + and - frequencies using exp(i*f*t), I use sines and cosines to 
% keep tc real.  If the input series is real, than this will 
% automatically use real-only computation (faster). However, for the analysis, 
% it's handy to get the + and - frequencies ('ap' and 'am'), and so 
% that's what we do afterwards.

% The basic code solves the matrix problem Ac=x+errors where the functions to
% use in the fit fill up the A matrix, which is of size (number points)x(number
% constituents). This can get very, very large for long time series, and
% for this the more complex block processing algorithm was added. It should
% give identical results (up to roundoff error)

if strcmp(lsq(1:3),'dir'),

  if secular(1:3)=='lin',
    tc=[ones(length(t),1) cos((2*pi)*t*fu') sin((2*pi)*t*fu') t*(2/dt/nobsu)];
  else
    tc=[ones(length(t),1) cos((2*pi)*t*fu') sin((2*pi)*t*fu') ];
  end;
  
  coef=tc(gd,:)\xin(gd);

  z0=coef(1);
  ap=(coef(2:(1+mu))-i*coef((2+mu):(1+2*mu)))/2;  % a+ amplitudes
  am=(coef(2:(1+mu))+i*coef((2+mu):(1+2*mu)))/2;  % a- amplitudes
  if secular(1:3)=='lin',
    dz0=coef(end);
  else
    dz0=0;
  end;    
  xout=tc*coef;  % This is the time series synthesized from the analysis

else  % More complicated code required for long time series when memory may be
      % a problem. Modified from code submitted by Derek Goring (NIWA Chrischurch)
      
      % Basically the normal equations are formed (rather than using Matlab's \
      % algorithm for least squares); this can be done by adding up subblocks
      % of data. Notice how the code is messier, and we have to recalculate everything
      % to get the original fit.

  nsub=5000;  % Block length - doesn't matter really but should be small enough to
              % get allocated quickly	      
  if secular(1:3)=='lin',
    lhs=zeros(mu*2+2,mu*2+2); rhs=zeros(mu*2+2,1);
    for j1=1:nsub:ngood
      j2=min(j1 + nsub - 1,ngood);
      E=[ones(j2-j1+1,1) cos((2*pi)*t(gd(j1:j2))*fu') sin((2*pi)*t(gd(j1:j2))*fu') t(gd(j1:j2))*(2/dt/nobsu)];
      rhs=rhs + E'*xin(gd(j1:j2));
      lhs=lhs + E'*E;
    end;
  else  
    lhs=zeros(mu*2+1,mu*2+1); rhs=zeros(mu*2+1,1);
    for j1=1:nsub:ngood
      j2=min(j1 + nsub - 1,ngood);
      E=[ones(j2-j1+1,1) cos((2*pi)*t(gd(j1:j2))*fu') sin((2*pi)*t(gd(j1:j2))*fu')];
      rhs=rhs + E'*xin(gd(j1:j2));
      lhs=lhs + E'*E;
    end;
  end;
    
  coef=lhs\rhs;
  
  z0=coef(1);
  ap=(coef(2:(1+mu))-i*coef((2+mu):(1+2*mu)))/2;  % a+ amplitudes
  am=(coef(2:(1+mu))+i*coef((2+mu):(1+2*mu)))/2;  % a- amplitudes
  if secular(1:3)=='lin',
    dz0=coef(end);
  else
    dz0=0;
  end; 
  
  xout=xin; % Copies over NaN   
  if secular(1:3)=='lin',
    for j1=1:nsub:nobs
      j2=min(j1 + nsub - 1,nobs);
      E=[ones(j2-j1+1,1) cos((2*pi)*t(j1:j2)*fu') sin((2*pi)*t(j1:j2)*fu') t(j1:j2)*(2/dt/nobsu)];
      xout(j1:j2)=E*coef;
    end;
  else  
    for j1=1:nsub:nobs
      j2=min(j1 + nsub - 1,nobs);
      E=[ones(j2-j1+1,1) cos((2*pi)*t(j1:j2)*fu') sin((2*pi)*t(j1:j2)*fu')];
      xout(j1:j2)=E*coef;
    end;
  end;

end;

   
 
%----------------------------------------------------------------------
% Check variance explained (but do this with the original fit).

xres=xin-xout; % and the residuals!

if isreal(xin),    % Real time series
  varx=cov(xin(gd));varxp=cov(xout(gd));varxr=cov(xres(gd));
  fprintf('   percent of var residual after lsqfit/var original: %5.2f %%\n',100*(varxr/varx));  
else               % Complex time series
  varx=cov(real(xin(gd)));varxp=cov(real(xout(gd)));varxr=cov(real(xres(gd)));
  fprintf('   percent of X var residual after lsqfit/var original: %5.2f %%\n',100*(varxr/varx));

  vary=cov(imag(xin(gd)));varyp=cov(imag(xout(gd)));varyr=cov(imag(xres(gd)));
  fprintf('   percent of Y var residual after lsqfit/var original: %5.2f %%\n',100*(varyr/vary));
end;


%---------- Correct for prefiltering-----------------------------------

corrfac=interp1(corr_fs,corr_fac,fu);
% To stop things blowing up!
corrfac(corrfac>100 | corrfac <.01 | isnan(corrfac))=1;

ap=ap.*corrfac;
am=am.*conj(corrfac);

%---------------Nodal Corrections-------------------------------------- 						   
% Generate nodal corrections and calculate phase relative to Greenwich. 						   
% Note that this is a slightly weird way to do the nodal corrections,							   
% but is 'traditional'.  The "right" way would be to change the basis							   
% functions used in the least-squares fit above.									   

if ~isempty(lat) & ~isempty(stime),   % Time and latitude								   

  % Get nodal corrections at midpoint time.										   
  [v,u,f]=t_vuf(ltype,centraltime,[ju;jinf],lat);									   

  vu=(v+u)*360; % total phase correction (degrees)									   
  nodcor=['Greenwich phase computed with nodal corrections applied to amplitude \n and phase relative to center time'];    
elseif ~isempty(stime),    % Time only  										   
  % Get nodal corrections at midpoint time										   
  [v,u,f]=t_vuf(ltype,centraltime,[ju;jinf]);										   
  vu=(v+u)*360; % total phase correction (degrees)									   
  nodcor=['Greenwich phase computed, no nodal corrections'];								   
else   % No time, no latitude												   
  vu=zeros(length(ju)+length(jinf),1);  										   
  f=ones(length(ju)+length(jinf),1);											   
   nodcor=['Phases at central time'];											   
end															   
fprintf(['   ',nodcor,'\n']);												   


%---------------Inference Corrections----------------------------------
% Once again, the "right" way to do this would be to change the basis
% functions.
ii=find(isfinite(jref));
if ii,
  fprintf('   Do inference corrections\n');
  snarg=nobsu*pi*(fi(ii)   -fu(jref(ii)) )*dt;
  scarg=sin(snarg)./snarg;
 
  if size(inf.amprat,2)==1,    % For real time series
    pearg=     2*pi*(vu(mu+ii)-vu(jref(ii))+inf.ph(ii))/360;
    pcfac=inf.amprat(ii).*f(mu+ii)./f(jref(ii)).*exp(i*pearg);
    pcorr=1+pcfac.*scarg;
    mcfac=conj(pcfac);
    mcorr=conj(pcorr);
  else                          % For complex time series
    pearg=     2*pi*(vu(mu+ii)-vu(jref(ii))+inf.ph(ii,1))/360;
    pcfac=inf.amprat(ii,1).*f(mu+ii)./f(jref(ii)).*exp(i*pearg);
    pcorr=1+pcfac.*scarg;
    mearg=    -2*pi*(vu(mu+ii)-vu(jref(ii))+inf.ph(ii,2))/360;
    mcfac=inf.amprat(ii,2).*f(mu+ii)./f(jref(ii)).*exp(i*mearg);
    mcorr=1+mcfac.*scarg;
  end;
    
  ap(jref(ii))=ap(jref(ii))./pcorr;   % Changes to existing constituents
  ap=[ap;ap(jref(ii)).*pcfac];        % Inferred constituents

  am(jref(ii))=am(jref(ii))./mcorr;
  am=[am;am(jref(ii)).*mcfac];

  fu=[fu;fi(ii)];
  nameu=[nameu;namei(ii,:)];
end;

% --------------Error Bar Calculations---------------------------------
%
% Error bar calcs involve two steps:
%      1) Estimate the uncertainties in the analyzed amplitude
%         for both + and - frequencies (i.e., in 'ap' and 'am').
%         A simple way of doing this is to take the variance of the
%         original time series and divide it into the amount appearing
%         in the bandwidth of the analysis (approximately 1/length).
%         A more sophisticated way is to assume "locally white"
%         noise in the vicinity of, e.g., the diurnal consistuents.
%         This takes into account slopes in the continuum spectrum.
%
%      2) Transform those uncertainties into ones suitable for ellipse
%         parameters (axis lengths, angles). This can be done 
%         analytically for large signal-to-noise ratios. However, the 
%         transformation is non-linear at lows SNR, say, less than 10
%         or so.
%

xr=fixgaps(xres); % Fill in "internal" NaNs with linearly interpolated
                  % values so we can fft things.
nreal=1;

if strmatch(errcalc(2:end),'boot'),
  fprintf('   Using nonlinear bootstrapped error estimates\n');
  
  % "noise" matrices are created with the right covariance structure
  % to add to the analyzed components to create 'nreal' REPLICATES. 
  % 

  nreal=300;             % Create noise matrices 
  [NP,NM]=noise_realizations(xr(isfinite(xr)),fu,dt,nreal,errcalc);
      
  % All replicates are then transformed (nonlinearly) into ellipse 
  % parameters.  The computed error bars are then based on the std
  % dev of the replicates.

  AP=ap(:,ones(1,nreal))+NP;        % Add to analysis (first column
  AM=am(:,ones(1,nreal))+NM;        % of NM,NP=0 so first column of
                                    % AP/M holds ap/m).
  epsp=angle(AP)*180/pi;            % Angle/magnitude form:
  epsm=angle(AM)*180/pi;
  ap=abs(AP);
  am=abs(AM);
elseif strmatch(errcalc,'linear'),
  fprintf('   Using linearized error estimates\n');
  %
  % Uncertainties in analyzed amplitudes are computed in different
  % spectral bands. Real and imaginary parts of the residual time series
  % are treated separately (no cross-covariance is assumed).
  %
  % Noise estimates are then determined from a linear analysis of errors,
  % assuming that everything is uncorrelated. This is OK for scalar time
  % series but can fail for vector time series if the noise is not 
  % isotropic.
  
  [ercx,eicx]=noise_stats(xr(isfinite(xr)),fu,dt);
  % Note - here we assume that the error in the cos and sin terms is 
  % equal, and equal to total power in the encompassing frequency bin. 
  % It seems like there should be a factor of 2 here somewhere but it 
  % only works this way! <shrug>
  [emaj,emin,einc,epha]=errell(ap+am,i*(ap-am),ercx,ercx,eicx,eicx);

  epsp=angle(ap)*180/pi;
  epsm=angle(am)*180/pi;
  ap=abs(ap);
  am=abs(am);
else
  error(['Unrecognized type of error analysis: ''' errcalc ''' specified!']);
end;

%-----Convert complex amplitudes to standard ellipse parameters--------

aap=ap./f(:,ones(1,nreal));	% Apply nodal corrections and
aam=am./f(:,ones(1,nreal));	% compute ellipse parameters.

fmaj=aap+aam;                   % major axis
fmin=aap-aam;                   % minor axis

gp=mod( vu(:,ones(1,nreal))-epsp ,360); % pos. Greenwich phase in deg.
gm=mod( vu(:,ones(1,nreal))+epsm ,360); % neg. Greenwich phase in deg.

finc= (epsp+epsm)/2;
finc(:,1)=mod( finc(:,1),180 ); % Ellipse inclination in degrees
				% (mod 180 to prevent ambiguity, i.e., 
				% we always ref. against northern 
				% semi-major axis.
	
finc=cluster(finc,180); 	% Cluster angles around the 'true' 
                                % angle to avoid 360 degree wraps.

pha=mod( gp+finc ,360); 	% Greenwich phase in degrees.

pha=cluster(pha,360);		% Cluster angles around the 'true' angle
				% to avoid 360 degree wraps.

%----------------Generate 95% CI---------------------------------------
%% For bootstrapped errors, we now compute limits of the distribution.
if strmatch(errcalc(2:end),'boot'),
     %% std dev-based estimates.
     % The 95% CI are computed from the sigmas
     % by a 1.96 fudge factor (infinite degrees of freedom).
     % emaj=1.96*std(fmaj,0,2);
     % emin=1.96*std(fmin,0,2);
     % einc=1.96*std(finc,0,2);
     % epha=1.96*std(pha ,0,2);
     %% Median-absolute-deviation (MAD) based estimates.
     % (possibly more stable?)
      emaj=median(abs(fmaj-median(fmaj,2)*ones(1,nreal)),2)/.6375*1.96;
      emin=median(abs(fmin-median(fmin,2)*ones(1,nreal)),2)/.6375*1.96;
      einc=median(abs(finc-median(finc,2)*ones(1,nreal)),2)/.6375*1.96;
      epha=median(abs( pha-median( pha,2)*ones(1,nreal)),2)/.6375*1.96;
else
   % In the linear analysis, the 95% CI are computed from the sigmas
   % by this fudge factor (infinite degrees of freedom).
   emaj=1.96*emaj;
   emin=1.96*emin;
   einc=1.96*einc;
   epha=1.96*epha;
end;
				  
if isreal(xin),
    tidecon=[fmaj(:,1),emaj,pha(:,1),epha];
else
    tidecon=[fmaj(:,1),emaj,fmin(:,1),emin, finc(:,1),einc,pha(:,1),epha];
end;

% Sort results by frequency (needed if anything has been inferred since 
% these are stuck at the end of the list by code above).
if any(isfinite(jref)),
 [fu,I]=sort(fu);
 nameu=nameu(I,:);
 tidecon=tidecon(I,:);
end;

snr=(tidecon(:,1)./tidecon(:,2)).^2;  % signal to noise ratio

%--------Generate a 'prediction' using significant constituents----------
xoutOLD=xout;
if synth>=0,
 if ~isempty(lat) & ~isempty(stime),
   fprintf('   Generating prediction with nodal corrections, SNR is %f\n',synth);
   xout=t_predic(stime+[0:nobs-1]*dt/24.0,nameu,fu,tidecon,'lat',lat,'synth',synth,'anal',ltype);
 elseif ~isempty(stime), 
   fprintf('   Generating prediction without nodal corrections, SNR is %f\n',synth);
   xout=t_predic(stime+[0:nobs-1]*dt/24.0,nameu,fu,tidecon,'synth',synth,'anal',ltype);
 else
   fprintf('   Generating prediction without nodal corrections, SNR is %f\n',synth);
   xout=t_predic(t/24.0,nameu,fu,tidecon,'synth',synth,'anal',ltype);
 end;
else
 fprintf('   Returning fitted prediction\n');
end;

%----------------------------------------------------------------------
% Check variance explained (but now do this with the synthesized fit).
xres=xin(:)-xout(:); % and the residuals!

%error;

if isreal(xin),    % Real time series
  varx=cov(xin(gd));varxp=cov(xout(gd));varxr=cov(xres(gd));
  fprintf('   percent of var residual after synthesis/var original: %5.2f %%\n',100*(varxr/varx));  
else               % Complex time series
  varx=cov(real(xin(gd)));varxp=cov(real(xout(gd)));varxr=cov(real(xres(gd)));
  fprintf('   percent of X var residual after synthesis/var original: %5.2f %%\n',100*(varxr/varx));

  vary=cov(imag(xin(gd)));varyp=cov(imag(xout(gd)));varyr=cov(imag(xres(gd)));
  fprintf('   percent of Y var residual after synthesis/var original: %5.2f %%\n',100*(varyr/vary));
end;


%-----------------Output results---------------------------------------

if fid>1,
 fprintf(fid,'\n%s\n',['file name: ',filen]);
elseif fid==1,
 fprintf(fid,'-----------------------------------\n');
end

if fid>0,
  fprintf(fid,'date: %s\n',date);
  fprintf(fid,'nobs = %d,  ngood = %d,  record length (days) = %.2f\n',nobs,ngood,length(xin)*dt/24);
  if ~isempty(stime); fprintf(fid,'%s\n',['start time: ',datestr(stime)]); end
  fprintf(fid,'rayleigh criterion = %.1f\n',ray);
  fprintf(fid,'%s\n',nodcor);
%  fprintf(fid,'\n     coefficients from least squares fit of x\n');
%  fprintf(fid,'\n tide    freq        |a+|       err_a+      |a-|       err_a-\n');
%  for k=1:length(fu);
%    if ap(k)>eap(k) | am(k)>eam(k), fprintf('*'); else fprintf(' '); end;
%    fprintf(fid,'%s  %8.5f  %9.4f  %9.4f  %9.4f  %9.4f\n',nameu(k,:),fu(k),ap(k),eap(k),am(k),eam(k));
%  end
  fprintf(fid,'\nx0= %.3g, x trend= %.3g\n',real(z0),real(dz0));
  fprintf(fid,['\nvar(x)= ',num2str(varx),'   var(xp)= ',num2str(varxp),'   var(xres)= ',num2str(varxr) '\n']);
  fprintf(fid,'percent var predicted/var original= %.1f %%\n',100*varxp/varx);

  if isreal(xin)
    fprintf(fid,'\n     tidal amplitude and phase with 95%% CI estimates\n');
    fprintf(fid,'\ntide   freq       amp     amp_err    pha    pha_err     snr\n');
    for k=1:length(fu);
      if snr(k)>synth, fprintf(fid,'*'); else fprintf(fid,' '); end;
      fprintf(fid,'%s %9.7f %9.4f %8.3f %8.2f %8.2f %8.2g\n',nameu(k,:),fu(k),tidecon(k,:),snr(k));
    end
  else
    fprintf(fid,'\ny0= %.3g, x trend= %.3g\n',imag(z0),imag(dz0));
    fprintf(fid,['\nvar(y)= ',num2str(vary),'    var(yp)= ',num2str(varyp),'  var(yres)= ',num2str(varyr) '\n']);
    fprintf(fid,'percent var predicted/var original= %.1f %%\n',100*varyp/vary);
    fprintf(fid,'\n%s\n',['ellipse parameters with 95%% CI estimates']);
    fprintf(fid,'\n%s\n',['tide   freq      major  emaj    minor   emin     inc    einc     pha    epha      snr']);
    for k=1:length(fu);
      if snr(k)>synth, fprintf(fid,'*'); else fprintf(fid,' '); end;
      fprintf(fid,'%s %9.7f %6.3f %7.3f %7.3f %6.2f %8.2f %6.2f %8.2f %6.2f %6.2g\n',...
	nameu(k,:),fu(k),tidecon(k,:),snr(k));
    end
    fprintf(fid,['\ntotal var= ',num2str(varx+vary),'   pred var= ',num2str(varxp+varyp) '\n']);
    fprintf(fid,'percent total var predicted/var original= %.1f %%\n\n',100*(varxp+varyp)/(varx+vary));
  end

  if fid~=1, st=fclose(fid); end
end;

xout=reshape(xout,inn,inm);
switch nargout,
  case {0,3,4}
  case {1}
   nameu = struct('name',nameu,'freq',fu,'tidecon',tidecon,'type',ltype);
  case {2}   
   nameu = struct('name',nameu,'freq',fu,'tidecon',tidecon,'type',ltype);
   fu=xout;
end;
   
%----------------------------------------------------------------------
function [nameu,fu,ju,namei,fi,jinf,jref]=constituents(minres,constit,...
                                     shallow,infname,infref,centraltime);
% [name,freq,kmpr]=constituents(minres,infname) loads tidal constituent
% table (containing 146 constituents), then picks out only the '
% resolvable' frequencies (i.e. those that are MINRES apart), base on 
% the comparisons in the third column of constituents.dat. Only 
% frequencies in the 'standard' set of 69 frequencies are actually used.
% Also return the indices of constituents to be inferred.

% If we have the mat-file, read it in, otherwise create it and read
% it in!

% R Pawlowicz 9/1/01 
% Version 1.0
%
%    19/1/02 - typo fixed (thanks to  Zhigang Xu)

% Compute frequencies from astronomical considerations.

if minres>1/(18.6*365.25*24),                       % Choose only resolveable pairs for short
  [const,sat,cshallow]=t_getconsts(centraltime);    % Time series   
  ju=find(const.df>=minres);
else                                                % Choose them all if > 18.6 years.
  [const,sat,cshallow]=t_get18consts(centraltime);
  ju=[2:length(const.freq)]';  % Skip Z0
  for ff=1:2,                  % loop twice to make sure of neightbouring pairs
    jck=find(diff(const.freq(ju))<minres);
    if (length(jck)>0)
       jrm=jck;
       jrm=jrm+(abs(const.doodsonamp(ju(jck+1)))<abs(const.doodsonamp(ju(jck))));
       disp(['  Warning! Following constituent pairs violate Rayleigh criterion']);
       for ick=1:length(jck);
	 disp(['     ',const.name(ju(jck(ick)),:),' vs ',const.name(ju(jck(ick)+1),:) ' - not using ',const.name(ju(jrm(ick)),:)]);
       end;
       ju(jrm)=[];
    end
  end;
end;
  
if ~isempty(constit),     % Selected if constituents are specified in input.
  ju=[];
  for k=1:size(constit,1),
   j1=strmatch(constit(k,:),const.name);
   if isempty(j1),
     disp(['Can''t recognize name ' constit(k,:) ' for forced search']);
   elseif j1==1,
     disp(['*************************************************************************']);
     disp(['Z0 specification ignored - for non-tidal offsets see ''secular'' property']);
     disp(['*************************************************************************']);
   else  
     ju=[ju;j1];
   end;
  end;
  [dum,II]=sort(const.freq(ju)); % sort in ascending order of frequency.
  ju=ju(II);
end;


disp(['   number of standard constituents used: ',int2str(length(ju))])

if ~isempty(shallow),          % Add explictly selected shallow water constituents.
 for k=1:size(shallow,1),
   j1=strmatch(shallow(k,:),const.name);
   if isempty(j1),
     disp(['Can''t recognize name ' shallow(k,:) ' for forced search']);
   else
     if isnan(const.ishallow(j1)),
       disp([shallow(k,:) ' Not a shallow-water constituent']);
     end;
     disp(['   Forced fit to ' shallow(k,:)]);
     ju=[ju;j1];
   end;
 end;
 
end;
      
nameu=const.name(ju,:);
fu=const.freq(ju);


% Check if neighboring chosen constituents violate Rayleigh criteria.
jck=find(diff(fu)<minres);
if (length(jck)>0)
   disp(['  Warning! Following constituent pairs violate Rayleigh criterion']);
   for ick=1:length(jck);
   disp(['     ',nameu(jck(ick),:),'  ',nameu(jck(ick)+1,:)]);
   end;
end

% For inference, add in list of components to be inferred.

fi=[];namei=[];jinf=[];jref=[];
if ~isempty(infname),
  fi=zeros(size(infname,1),1);
  namei=zeros(size(infname,1),4);
  jinf=zeros(size(infname,1),1)+NaN;
  jref=zeros(size(infname,1),1)+NaN;

  for k=1:size(infname,1),
   j1=strmatch(infname(k,:),const.name);
   if isempty(j1),
     disp(['Can''t recognize name' infname(k,:) ' for inference']);
   else
    jinf(k)=j1;
    fi(k)=const.freq(j1);
    namei(k,:)=const.name(j1,:);
    j1=strmatch(infref(k,:),nameu);
    if isempty(j1),
      disp(['Can''t recognize name ' infref(k,:) ' for as a reference for inference']);
    else
      jref(k)=j1;
      fprintf(['   Inference of ' namei(k,:) ' using ' nameu(j1,:) '\n']);
    end;
   end;
  end;    
  jinf(isnan(jref))=NaN;
end;

%----------------------------------------------------------------------
function y=fixgaps(x);
% FIXGAPS: Linearly interpolates gaps in a time series
% YOUT=FIXGAPS(YIN) linearly interpolates over NaN in the input time 
% series (may be complex), but ignores trailing and leading NaNs.

% R. Pawlowicz 11/6/99
% Version 1.0

y=x;

bd=isnan(x);
gd=find(~bd);

bd([1:(min(gd)-1) (max(gd)+1):end])=0;


y(bd)=interp1(gd,x(gd),find(bd)); 


%----------------------------------------------------------------------
function ain=cluster(ain,clusang);
% CLUSTER: Clusters angles in rows around the angles in the first 
% column. CLUSANG is the allowable ambiguity (usually 360 degrees but
% sometimes 180).

ii=(ain-ain(:,ones(1,size(ain,2))))>clusang/2;
ain(ii)=ain(ii)-clusang;
ii=(ain-ain(:,ones(1,size(ain,2))))<-clusang/2;
ain(ii)=ain(ii)+clusang;


%----------------------------------------------------------------------
function [NP,NM]=noise_realizations(xres,fu,dt,nreal,errcalc);
% NOISE_REALIZATIONS: Generates matrices of noise (with correct
% cross-correlation structure) for bootstrap analysis.
%

% R. Pawlowicz 11/10/00
% Version 1.0

if strmatch(errcalc,'cboot'),
  [fband,Pxrave,Pxiave,Pxcave]=residual_spectrum(xres,fu,dt);
  
  Pxcave=zeros(size(Pxcave));  %% For comparison with other technique!
  %fprintf('**** Assuming no covariance between u and v errors!*******\n');

elseif strmatch(errcalc,'wboot'),
  fband=[0 .5];
  nx=length(xres);
  A=cov(real(xres),imag(xres))/nx;
  Pxrave=A(1,1);Pxiave=A(2,2);Pxcave=A(1,2);
else
  error(['Unrecognized type of bootstap analysis specified: ''' errcalc '''']);
end;
  
nfband=size(fband,1);

Mat=zeros(4,4,nfband);
for k=1:nfband,

  % The B matrix represents the covariance matrix for the vector
  % [Re{ap} Im{ap} Re{am} Im{am}]' where Re{} and Im{} are real and
  % imaginary parts, and ap/m represent the complex constituent 
  % amplitudes for positive and negative frequencies when the input
  % is bivariate white noise. For a flat residual spectrum this works 
  % fine.
 
  % This is adapted here for "locally white" conditions, but I'm still
  % not sure how to handle a complex sxy, so this is set to zero
  % right now.
  
  p=(Pxrave(k)+Pxiave(k))/2;
  d=(Pxrave(k)-Pxiave(k))/2;
  sxy=Pxcave(k);
  
  B=[p    0   d   sxy;
     0    p  sxy  -d;
     d   sxy  p    0
     sxy -d   0    p];

  % Compute the transformation matrix that takes uncorrelated white 
  % noise and makes noise with the same statistical structure as the 
  % Fourier transformed noise.
  [V,D]=eig(B);
  Mat(:,:,k)=V*diag(sqrt(diag(D)));
end;

% Generate realizations for the different analyzed constituents.

N=zeros(4,nreal);
NM=zeros(length(fu),nreal);
NP=NM;
for k=1:length(fu);
  l=find(fu(k)>fband(:,1) & fu(k)<fband(:,2));
  N=[zeros(4,1),Mat(:,:,l)*randn(4,nreal-1)];
  NP(k,:)=N(1,:)+i*N(2,:);
  NM(k,:)=N(3,:)+i*N(4,:);
end;

%----------------------------------------------------------------------
function [ercx,eicx]=noise_stats(xres,fu,dt);
% NOISE_STATS: Computes statistics of residual energy for all 
% constituents (ignoring any cross-correlations between real and
% imaginary parts).

% S. Lentz  10/28/99
% R. Pawlowicz 11/1/00
% Version 1.0

[fband,Pxrave,Pxiave,Pxcave]=residual_spectrum(xres,fu,dt);
nfband=size(fband,1);
mu=length(fu);

% Get the statistics for each component.
ercx=zeros(mu,1);
eicx=zeros(mu,1);
for k1=1:nfband;
   k=find(fu>=fband(k1,1) & fu<=fband(k1,2));
   ercx(k)=sqrt(Pxrave(k1));
   eicx(k)=sqrt(Pxiave(k1));
end

%----------------------------------------------------------------------
function [fband,Pxrave,Pxiave,Pxcave]=residual_spectrum(xres,fu,dt)
% RESIDUAL_SPECTRUM: Computes statistics from an input spectrum over
% a number of bands, returning the band limits and the estimates for
% power spectra for real and imaginary parts and the cross-spectrum.          
%
% Mean values of the noise spectrum are computed for the following 
% 8 frequency bands defined by their center frequency and band width:
% M0 +.1 cpd; M1 +-.2 cpd; M2 +-.2 cpd; M3 +-.2 cpd; M4 +-.2 cpd; 
% M5 +-.2 cpd; M6 +-.21 cpd; M7 (.26-.29 cpd); and M8 (.30-.50 cpd). 

% S. Lentz  10/28/99
% R. Pawlowicz 11/1/00
% Version 1.0

% Define frequency bands for spectral averaging.
fband =[.00010 .00417;
        .03192 .04859;
        .07218 .08884;
        .11243 .12910;
        .15269 .16936;
        .19295 .20961;
        .23320 .25100;
        .26000 .29000;
        .30000 .50000];

% If we have a sampling interval> 1 hour, we might have to get
% rid of some bins.
%fband(fband(:,1)>1/(2*dt),:)=[];

nfband=size(fband,1);
nx=length(xres);

% Spectral estimate (takes real time series only).


% Matlab has changed their spectral estimator functions
% To match the old code, I have to divide by 2*dt. This is because
% 
%  PSD*dt  is two-sided spectrum in units of power per hertz.
%
%  PWELCH is the one-sided spectrum in power per hertz
%
%  So PWELCH/2 = PSD*dt


%[Pxr,fx]=psd(real(xres),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme. If you have an error here you are probably missing this toolbox
%[Pxi,fx]=psd(imag(xres),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
%[Pxc,fx]=csd(real(xres),imag(xres),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.


[Pxr,fx]=pwelch(real(xres),hanning(nx),ceil(nx/2),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme. If you have an error here you are probably missing this toolbox
Pxr=Pxr/2/dt;
[Pxi,fx]=pwelch(imag(xres),hanning(nx),ceil(nx/2),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
Pxi=Pxi/2/dt;
[Pxc,fx]=cpsd(real(xres),imag(xres),[],[],nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
Pxc=Pxc/2/dt;

df=fx(3)-fx(2);
Pxr(round(fu./df)+1)=NaN ; % Sets Px=NaN in bins close to analyzed frequencies
Pxi(round(fu./df)+1)=NaN ; % (to prevent leakage problems?).
Pxc(round(fu./df)+1)=NaN ; 

Pxrave=zeros(nfband,1);
Pxiave=zeros(nfband,1);
Pxcave=zeros(nfband,1);
% Loop downwards in frequency through bands (cures short time series
% problem with no data in lowest band).
%
% Divide by nx to get power per frequency bin, and multiply by 2
% to account for positive and negative frequencies.
%
for k=nfband:-1:1,
   jband=find(fx>=fband(k,1) & fx<=fband(k,2) & isfinite(Pxr));
   if any(jband),
     Pxrave(k)=mean(Pxr(jband))*2/nx;
     Pxiave(k)=mean(Pxi(jband))*2/nx;
     Pxcave(k)=mean(Pxc(jband))*2/nx;
   elseif k<nfband,
     Pxrave(k)=Pxrave(k+1);   % Low frequency bin might not have any points...
     Pxiave(k)=Pxiave(k+1);   
     Pxcave(k)=Pxcave(k+1);   
   end;
end
 

%----------------------------------------------------------------------
function [emaj,emin,einc,epha]=errell(cxi,sxi,ercx,ersx,ercy,ersy)
% [emaj,emin,einc,epha]=errell(cx,sx,cy,sy,ercx,ersx,ercy,ersy) computes
% the uncertainities in the ellipse parameters based on the 
% uncertainities in the least square fit cos,sin coefficients.
%
%  INPUT:  cx,sx=cos,sin coefficients for x 
%          cy,sy=cos,sin coefficients for y
%          ercx,ersx=errors in x cos,sin coefficients
%          ercy,ersy=errors in y cos,sin coefficients
%          
%  OUTPUT: emaj=major axis error
%          emin=minor axis error
%          einc=inclination error (deg)
%          epha=pha error (deg)

% based on linear error propagation, with errors in the coefficients 
% cx,sx,cy,sy uncorrelated. 

% B. Beardsley  1/15/99; 1/20/99
% Version 1.0

r2d=180./pi;
cx=real(cxi(:));sx=real(sxi(:));cy=imag(cxi(:));sy=imag(sxi(:));
ercx=ercx(:);ersx=ersx(:);ercy=ercy(:);ersy=ersy(:);

rp=.5.*sqrt((cx+sy).^2+(cy-sx).^2);
rm=.5.*sqrt((cx-sy).^2+(cy+sx).^2);
ercx2=ercx.^2;ersx2=ersx.^2;
ercy2=ercy.^2;ersy2=ersy.^2;

% major axis error
ex=(cx+sy)./rp;
fx=(cx-sy)./rm;
gx=(sx-cy)./rp;
hx=(sx+cy)./rm;
dcx2=(.25.*(ex+fx)).^2;
dsx2=(.25.*(gx+hx)).^2;
dcy2=(.25.*(hx-gx)).^2;
dsy2=(.25.*(ex-fx)).^2;
emaj=sqrt(dcx2.*ercx2+dsx2.*ersx2+dcy2.*ercy2+dsy2.*ersy2);

% minor axis error
dcx2=(.25.*(ex-fx)).^2;
dsx2=(.25.*(gx-hx)).^2;
dcy2=(.25.*(hx+gx)).^2;
dsy2=(.25.*(ex+fx)).^2;
emin=sqrt(dcx2.*ercx2+dsx2.*ersx2+dcy2.*ercy2+dsy2.*ersy2);

% inclination error
rn=2.*(cx.*cy+sx.*sy);
rd=cx.^2+sx.^2-(cy.^2+sy.^2);
den=rn.^2+rd.^2;
dcx2=((rd.*cy-rn.*cx)./den).^2;
dsx2=((rd.*sy-rn.*sx)./den).^2;
dcy2=((rd.*cx+rn.*cy)./den).^2;
dsy2=((rd.*sx+rn.*sy)./den).^2;
einc=r2d.*sqrt(dcx2.*ercx2+dsx2.*ersx2+dcy2.*ercy2+dsy2.*ersy2);

% phase error
rn=2.*(cx.*sx+cy.*sy);
rd=cx.^2-sx.^2+cy.^2-sy.^2;
den=rn.^2+rd.^2;
dcx2=((rd.*sx-rn.*cx)./den).^2;
dsx2=((rd.*cx+rn.*sx)./den).^2;
dcy2=((rd.*sy-rn.*cy)./den).^2;
dsy2=((rd.*cy+rn.*sy)./den).^2;
epha=r2d.*sqrt(dcx2.*ercx2+dsx2.*ersx2+dcy2.*ercy2+dsy2.*ersy2);





