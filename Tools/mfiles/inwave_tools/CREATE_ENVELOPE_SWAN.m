%close all
%clear all

% add wafo tooldbox

%cd 'C:\Program Files\MATLAB\R2007b\toolbox\wafo-2.1.1'
%initwafo

%cd 'K:\ISABEL_INWAVE\MATLAB_SCRIPTS'

% ANALYSIS PERIOD

%dateini=datenum(2003,09,16,12,0,0);
%datefin=datenum(2003,09,19,0,0,0);

%point=11;
%filename=strcat('K:\ISABEL_INWAVE\OUTPUTS\run_8\point',num2str(point),'.spc2d');
%filename='point1.spc2d';

filename='D:\Sandy_InWave\offbreach_2012103100.spc2d';
[FA1,DIR1,S1]=read_swan(filename);



%count=0;
%for ii=38:38+23

    figure
    SSS(:,:)=S1(:,:);
    pcolor(FA1,DIR1,SSS)
    shading flat
%   axis([0 0.3 100 180])
    caxis([0 10])
    ylabel('Direction, degrees')
    xlabel('Frequency (s^{-1})')
    colorbar 'vertical'
    figure1=(strcat('espectra_',num2str(1)));
%    print('-dpng','-r200',figure1)

    DIR1(DIR1<0)=360+DIR1(DIR1<0);
%   DIR1=DIR1;

    swan_f=FA1(:,1);
    swan_dir=DIR1(1,:)*pi/180;
    SSS=SSS./180;  %units J/m2/Hz/rad

    dum_dir=find(SSS==max(SSS(:)));
    [fpos,dirpos]=ind2ij(SSS,dum_dir);

    wave_dir=swan_dir(dirpos).*180/pi;
    TA=1./swan_f(fpos);

    df_swan=swan_f(2:end)-swan_f(1:end-1);
    ddir_swan=abs((swan_dir(2:end)-swan_dir(1:end-1)));
    SF=sum(SSS,2).*ddir_swan(1,1);
    m0 = trapz(swan_f,SF);
    Hs1=4.004*sqrt(m0);


    % INTERPOLATE TO A MORE RESOLUTION SPECTRA

    nf=100;
    ndir=60;
    Df=(0.5)/(nf);
    Ddir=(2*pi)/(ndir);

    % FREQUENCY LIMITS
    fi=Df/2;
    ff=0.5-Df/2;

    % DISCRETIZATION IN ENERGY AND DIRECTIONS
    [f,theta]=meshgrid(linspace(fi,ff,nf),linspace(Ddir,2*pi,ndir));

    % COMPUTATION OF THE SPECTRAL ENERGY

    Sp(:,:)=griddata(FA1,DIR1.*pi/180,SSS,f,theta);
    Sp(isnan(Sp)==1)=0;
    SW=sum(Sp*Ddir);

    % PLOT SPECTRA IN CARTESIAN SYSTEM

    figure(1)
    plot(f,SW,'linewidth',2)
    hold on
    plot(swan_f,SF,'r','linewidth',2)

    xlabel('Frequency (1/s)')
%   print('-dpng','-r200','spectral_ver')


    % COMPUTATION AF THE AMPLITUDE OF EACH FREQUENCY COMPONENT
    clear w
    w(:,1)=2*pi*f(1,:);
    Dw=2*pi*Df(1,:);
    aj1=zeros(nf);
    comp1=find(SW(1,:)>0);
    aj1(comp1)=sqrt(2*SW(1,comp1)*Df);
    figure5=('envelope_spectra_case');

    % RANDOM PHASE FOR EACH COMPONENT
    eps=(2*pi*rand(length(comp1))-pi);

    ind=0;
    time=[0:1:3600];
    eta1=zeros(length(time),1);

    %para 1000 s cada 0.05 s
    count=0;
    for t=0:1:3600
        count=count+1;
        ind=ind+1;
        for j=1:length(comp1)
            eta1(ind,1)=eta1(ind,1)+aj1(comp1(j)).*cos(-w(comp1(j)).*t+eps(comp1(j)));
        end
    end

    figure2=strcat('signal_case',num2str(1));
    figure2_2=strcat('signal_case_short',num2str(1));
    figure3=strcat('fft_envelope_case',num2str(1));
    figure4=strcat('reconstructed_spectra_case',num2str(1));


    % graphic representation of the generated serie

    figure (2)
    title('Free surface elevation')
    plot(time,eta1,'LineWidth',2)
    hold on
    xlabel('Time, sec')
    ylabel('\eta, m')

%    xx=[time; eta1'];
%    S=dat2spec2(xx',2000);
%    Sf = ttspec(S,'f');

%    Df=(Sf.f(2,1)-Sf.f(1,1));
%    Hs_p=4.004*sqrt(sum(Sf.S).*Df);

    % OBTAINING THE ENVELOPE

    figure(2)
    x_hil=hilbert(eta1);
    amp_hil=sqrt(real(x_hil).^2+imag(x_hil).^2);
    plot(time,amp_hil,'g','LineWidth',2)
    axis([0 1000 -5 5])
    legend('Sea surface elevation','Envelope')
    hold off
%    print('-dpng','-r200',figure2)
%    axis([300 600 -5 5])
%    print('-dpng','-r200',figure2_2)

    t=time;
    A=amp_hil;
    df=1/(t(2)-t(1));
    a=fft(A);
    a=a';
    f1=[linspace(0,df/2,length(a)/2+1)];
    a=2*abs(a(1:length(f1)))/length(t);
    a1=2*abs(a(1:length(f1)))/length(t);

    figure(3)
    plot(1./f1,a,'.-')
    axis([0 500 0 1])
    xlabel('Period, s')
    title('ENVELOPE FFT')
    ylabel('Amplitude of the envelope,m')
%    print(gcf, '-djpeg', '-zbuffer', figure3);



%    xx2=[time', amp_hil-mean(amp_hil)];
%    S=dat2spec2(xx2,2000)
%    Sf_env = ttspec(S,'f');

%    figure(5)
%    hold on
%    ii1=ceil((ii-37)/2)
%    subplot(4,4,ii1)
%    if mod(ii-37,2)==0
%        plot(1./Sf_env.f,Sf_env.S,'LineWidth',2)
%    else
%        plot(1./Sf_env(ii-37).f,Sf_env(ii-37).S,'r','LineWidth',2)
%    end
%    axis([0 500 0 40])
%    if mod(ii-37,2)==0
%        legend(strcat('Case',num2str(ii-1-37)),strcat('Case',num2str(ii-37)))
%    end

%    xlabel('Period, (s)')
%    ylabel('Spectral energy density')
%    hold on



%    figure(5)
%    title('ENVELOPE SPECTRA')
%    set(gcf,'PaperPosition',[0.1 0.1 20 20]);
%    print(gcf, '-djpeg', '-zbuffer', figure5);


if(0)
  weigths=zeros(24,24*3600);
  eta=zeros(24,24*3600);

  for ii=1:24
    weigths(ii,(ii-1)*3600+1:(ii-1)*3600+3600)=1;
    if ii>1
        weigths(ii,(ii-1)*3600+1:(ii-1)*3600+300)=(0:1./(3900-3601):1);     end
    if ii<24
        weigths(ii,(ii)*3600+1:(ii)*3600+300)=(1:-1./(3900-3601):0);
    end
    if ii<24
    eta(ii,(ii-1)*3600+1:(ii-1)*3600+3900)=amp_hil(ii,1:end-1).*weigths(ii,(ii-1)*3600+1:(ii-1)*3600+3900);
    else
    eta(ii,(ii-1)*3600+1:(ii-1)*3600+3600)=amp_hil(ii,1:end-301).*weigths(ii,(ii-1)*3600+1:(ii-1)*3600+3600);
    end        
  end
end

eta=amp_hil;
eta_tot=eta;
%eta_tot=sum(eta);
%
%
%save realizations amp_hil time Sf_env a1 f1 TA wave_dir Hs1
save realizations amp_hil time a1 f1 Hs1 eta_tot
%
%
%