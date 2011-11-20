clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  THIS SCRIPT VERIFIES THE PROCEDURES IN INWAVE %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%    TO REPRODUCE WAVE CURRENT INTERACTION     %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) DEFINE TIME VECTOR, SIMULATION DOMAIN, WATER DEPTH and CURRENTS.

dt=600;
dx=1000;
dy=1000;

t=[0:dt:24*3600];

[X,Y]=meshgrid(0:dx:50*dx,0:dy:60*dy);

[nm,nn]=size(X);

nt=length(t);

h0=10*ones(size(X));

%u0=zeros(size(X));
u0=-exp( -((X-25000)./5000).^2-((60*1000-Y)./20000).^2);
pcolor(u0)
shading flat
colorbar 'vertical'
u=zeros(nt,nm,nn);
h=zeros(nt,nm,nn);

for tt=1:nt
    for jj=1:nm
        for ii=1:nn
            u(tt,jj,ii)=0.1*u0(jj,ii)*sin(2*pi/(12*3600).*t(1,tt));
            h(tt,jj,ii)=h0(jj,ii)+sin(2*pi/(12*3600).*t(1,tt));
        end
    end
end

% 2) DEFINE THE ABSOLUTE PERIOD

Ta=ones(size(h0)).*10;

% 3) DEFINE THE RELATIVE PERIOD

Ta;


wr=zeros(size(u));
wa=zeros(size(u));
Tr=zeros(size(u));
k=zeros(size(u));

Tr(1,:,:)=Ta;
wr(1,:,:)=(2*pi)./Tr(1,:,:);


wa=zeros(size(u));

clear L error

for tt=1:nt-1

display(tt) 
  
  for jj=1:nm
    for ii=1:nn
      error=1000;
      L0(jj,ii)=(9.81*Tr(tt,jj,ii)^2)./(2*pi);
      L1(jj,ii)=L0(jj,ii)*tanh(2*pi*h(tt,jj,ii)/L0(jj,ii));
      while(error>0.1)
        L(jj,ii)=L0(jj,ii)*tanh(2*pi*h(tt,jj,ii)./L1(jj,ii));
        error=abs(L(jj,ii)-L1(jj,ii))/L(jj,ii)*100;
        L1(jj,ii)=L(jj,ii);
      end
      k(tt,jj,ii)=(2*pi/L(jj,ii)); 
    end
  end


% 6) COMPUTE THE WAVE ABSOLUTE FREQUENCY
    
wa(tt,:,:)=wr(tt,:,:)+k(tt,:,:).*u(tt,:,:) ;


% 7) COMPUTE THE CHANGE IN THE WAVE NUMBER

for jj=2:nm-1
  for ii=1:nn
    k(tt+1,jj,ii)=k(tt,jj,ii)+(wa(tt,jj+1,ii)-wa(tt,jj-1,ii))./dy*dt;
  end
end

k(tt+1,1,:)=k(tt+1,2,:);
k(tt+1,nm,:)=k(tt+1,nm-1,ii);

% 8) COMPUTE THE RELATIVE FREQUENCY

wr(tt+1,:,:)=((9.81*k(tt+1,:,:)).*tanh(k(tt+1,:,:).*h(tt+1,:,:))).^0.5; 
Tr(tt+1,:,:)=2*pi/wr(tt+1,:,:);


    figure(1)
    
    subplot(2,2,1)
    pcolor(squeeze(h(tt,:,:)))
    shading flat
    colorbar 'vertical'
    title('h')

    subplot(2,2,2)
    pcolor(squeeze(u(tt,:,:)))
    shading flat
    colorbar 'vertical'
    caxis([-0.2 0.2])
    title('u')
    
    subplot(2,2,3)
    pcolor((2*pi)./squeeze(k(tt,:,:)))
    shading flat
    colorbar 'vertical'
    title('L')
    
    subplot(2,2,4)
    pcolor((2*pi)./squeeze(wr(tt,:,:)))
    shading flat
    colorbar 'vertical'
    caxis([9.5 10.5])
    title('wr')
    pause
  
end


figure
plot(2*pi./squeeze(wr(:,31,26)))
hold on
plot(2*pi./squeeze(wa(:,31,26)),'r')
    
