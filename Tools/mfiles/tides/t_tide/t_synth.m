function [sm,lm,tcon]=t_synth(varargin)
% T_SYNTH Monte-Carlo test of the error estimation using synthetic data
% A single test of the harmonic analysis involves:
%
%   1) Generation of a "pure" tidal signal with known constituents.
%   2) Contamination with noise (according to some statistical model).
%   3) The harmonic analysis resulting in constituent estimates.
%   4) Generation of confidence intervals for the estimates.
%
% T_SYNTH runs multiple realizations of this sequence. In each 
% realization, the added noise is different (but its statistical model
% remains the same).  Statistically, one would hope that the "true" 
% constituent components were somewhere within the 95% CI 95% of the 
% time. Almost equivalently, the width of the 95% CI should match the
% 95%-width of the histogram of estimates. The resulting plot shows 
% the histogram of estimates, their 95% width, and the width of 95% 
% CI for all realizations.
%
% A variety of input parameters can be specified:
%
%       'freqs':   Frequencies to use (default {'M2','K1','S2'} ).
%
%       'tidecon': Constituents in run for specified frequencies in
%                  the form of a matrix with one row per constituent.
%                  Each row is of the form 
%                   [semi-major_axis_length, semi-minor_axis_length,
%                    ellipse_inclination (deg), Greenwich_phase (deg)]'
%                  If only one row is given, the same parameters are
%                  used for all constituents.
%                  default: [1 .1 45 60]
%             
%                  If you want to test things for a "real" time series,
%                  (e.g., an elevation series), set the two middle
%                  parameters to 0 - e.g., [1 0 0 60].
%                            
%       'time':    Time axis (hours) default [0:24*60] 
%
%       'nrun':    Number of simulations (default 100);
%
%       'error':   Formula for errors to be added to simulation (to be
%                  used in an EVAL statement). If the size of the 
%                  synthesized data matrix is needed, replace with SY. 
%                  Examples (of correlated bivariate gaussian white 
%                  noise) are:
%
%                    '.5*randn(SY)'
%                    '(1*randn(SY)*1+1e-6*i*randn(SY))*exp(i*pi/4)'
%               
%                  To get coloured noise, use the implicit function
%                  'colrand': e.g. to get noise with a wavenumber slope
%                  of -1.1
%                   
%                    '.5*colrand(SY,-1.1)';
%
%       'boota':   Bootstrap analysis. Either
%                   'c': Assume coloured uncorrelated noise.
%                   'w': Assume white bivariate noise.
%
%   Outputs: While running different realizations, text is output to
%   the console. Upon completion, a figure is drawn, in which coloured
%   histograms (a different colour for each constituent) of the 
%   stimated values of all realizations are drawn. + and - 2.5% 
%   percentiles are indicated with thick dashed bars (these limits are
%   taken as 1.96 * the standard deviation). Uncertainties for each 
%   realization are shown by thin solid lines; these should lie on top
%   of the thicker dashed curves if the error analysis is correct.  The 
%   upper set of plots uses the bootstrapped confidence intervals, and 
%   the lower set a 'linear' analysis.
%   


% R. Pawlowicz 6/5/00
%             11/2/00 - Added linear analysis.
% Version 1.0

const=t_getconsts;

% Defaults

freqs={'M2','K1','S2'};
%   May Min Inc Gphase
tcon=[1 .1 45 60];
t=[0:24*60];
nrun=100;
errstr='.5*randn(SY)+.5*i*randn(SY)';
boota='c';

while length(varargin)>0,
  if isstr(varargin{1}),
    switch lower(varargin{1}(1:3)),
      case 'fre',
        freqs=varargin{2};
      case 'tid',
        tcon=varargin{2};
      case 'tim',
        t=varargin{2};
      case 'nru',
        nrun=varargin{2};
      case 'err',
        errstr=varargin{2};
      case 'boo',
        boota=varargin{2};
      otherwise,
        error(['Can''t understand property:' varargin{1}]);
      end;     
    varargin([1 2])=[];
  end;
end;


t=t-mean(t);

nrunreq=length(freqs);
freq=zeros(nrunreq,1);
for k=1:nrunreq,
  freq(k)=const.freq(strmatch(freqs(k),const.name));
end;

if size(tcon,1)==1, tcon=tcon(ones(nrunreq,1),:); end;
tscmplx=1;
if ~(any(tcon(:,2)) | any(tcon(:,3))), tscmplx=0;  end;

ap=(tcon(:,1)+tcon(:,2))/2;
am=(tcon(:,1)-tcon(:,2))/2;
em=(tcon(:,3)+tcon(:,4))*pi/180;
ep=(tcon(:,3)-tcon(:,4))*pi/180;


sm=zeros(nrunreq,8,nrun);
lm=zeros(nrunreq,8,nrun);
for k=1:nrun,

  y=sum(((ap.*exp(i*ep))*ones(1,length(t))).*exp(i*2*pi*freq*t)+...
	((am.*exp(i*em))*ones(1,length(t))).*exp(-i*2*pi*freq*t));
  SY=size(y);

  eval(['y=y+' errstr ';']);


  fprintf('%d/%d\n',k,nrun);

  [nameu,fu,tidecon,xout]=t_tide(y,'interval',t(2)-t(1),'output','none',...
                                 'rayleigh',freqs,'error',[boota 'boot']);
  [nameu2,fu2,tidecon2,xout2]=t_tide(y,'interval',t(2)-t(1),'output','none',...
                                 'rayleigh',freqs,'error','linear');

  I=zeros(nrunreq,1);
  for l=1:nrunreq,
    I(l)=strmatch(freqs(l),nameu);
  end;

  if tscmplx,
    sm(:,:,k)=tidecon(I,:);
    lm(:,:,k)=tidecon2(I,:);
  else
    sm(:,[1 2 7 8],k)=tidecon(I,:);
    lm(:,[1 2 7 8],k)=tidecon2(I,:);
  end;
end;
for k=1:size(sm,1),
  ii=sm(k,7,:)-tcon(k,4) > 180;
  sm(k,7,ii)=sm(k,7,ii)-360;
  ii=sm(k,7,:)-tcon(k,4) < -180;
  sm(k,7,ii)=sm(k,7,ii)+360;
  ii=lm(k,7,:)-tcon(k,4) > 180;
  lm(k,7,ii)=lm(k,7,ii)-360;
  ii=lm(k,7,:)-tcon(k,4) < -180;
  lm(k,7,ii)=lm(k,7,ii)+360;
end;

SM=sm; lbl='boot';
if boota=='c', lbl='Coloured Boot'; 
else           lbl='White Boot'; end;

for e=[0 4],

  if e==4, SM=lm; lbl='linear'; end;

  subplot(2,4,1+e);

  hist(squeeze(SM(:,1,:))')
  SS=1.96*std(squeeze(SM(:,1,:))')';
  MM=mean(squeeze(SM(:,1,:))')';
  line(([MM MM]+[SS SS])',[0 nrun/5],'linewidth',3,'linest','--');
  line(([MM MM]-[SS SS])',[0 nrun/5],'linewidth',3,'linest','--');
  line((tcon(:,ones(1,nrun))+squeeze(SM(:,2,:)))',[0:nrun-1]/5);
  line((tcon(:,ones(1,nrun))-squeeze(SM(:,2,:)))',[0:nrun-1]/5);
  title(['Umajor (' lbl ')']);

  subplot(2,4,2+e);

  hist(squeeze(SM(:,3,:))')
  SS=1.96*std(squeeze(SM(:,3,:))')';
  MM=mean(squeeze(SM(:,3,:))')';
  line(([MM MM]+[SS SS])',[0 nrun/5],'linewidth',3,'linest','--');
  line(([MM MM]-[SS SS])',[0 nrun/5],'linewidth',3,'linest','--');
  line((tcon(:,2*ones(1,nrun))+squeeze(SM(:,4,:)))',[0:nrun-1]/5);
  line((tcon(:,2*ones(1,nrun))-squeeze(SM(:,4,:)))',[0:nrun-1]/5);
  title(['Uminor (' lbl ')']);

  subplot(2,4,3+e);

  hist(squeeze(SM(:,5,:))')
  SS=1.96*std(squeeze(SM(:,5,:))')';
  MM=mean(squeeze(SM(:,5,:))')';
  line(([MM MM]+[SS SS])',[0 nrun/5],'linewidth',3,'linest','--');
  line(([MM MM]-[SS SS])',[0 nrun/5],'linewidth',3,'linest','--');
  line((tcon(:,3*ones(1,nrun))+squeeze(SM(:,6,:)))',[0:nrun-1]/5);
  line((tcon(:,3*ones(1,nrun))-squeeze(SM(:,6,:)))',[0:nrun-1]/5);
  title(['Inclination (' lbl ')']);
  xlabel('              --: Actual 95% CI     -: Estimated 95% CI','fontweight','bold');

  subplot(2,4,4+e);

  hist(squeeze(SM(:,7,:))')
  SS=1.96*std(squeeze(SM(:,7,:))')';
  MM=mean(squeeze(SM(:,7,:))')';
  line(([MM MM]+[SS SS])',[0 nrun/5],'linewidth',3,'linest','--');
  line(([MM MM]-[SS SS])',[0 nrun/5],'linewidth',3,'linest','--');
  line((tcon(:,4*ones(1,nrun))+squeeze(SM(:,8,:)))',[0:nrun-1]/5);
  line((tcon(:,4*ones(1,nrun))-squeeze(SM(:,8,:)))',[0:nrun-1]/5);
  title(['G-phase (' lbl ')']);

end;

if nargout==0,
 clear sm lm tcon
end;


%----------------------------------------------------------------------
function A=colrand(SY,slp)
% COLRAND generates coloured noise
% Version 1.0

lSY=max(SY);
f=[1 1:floor(lSY/2) -ceil(lSY/2)+1:-1]/(lSY);
x=randn(SY);
x2=(abs(f).^slp).*fft(x);
A=real(ifft(x2));
A=A/std(A);

%%plot([1:max(SY)]',[x' A']);





