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
            u(tt,jj,ii)=u0(jj,ii)*sin(2*pi/(12*3600).*t(1,tt));
            h(tt,jj,ii)=h0(jj,ii);  %+sin(2*pi/(12*3600).*t(1,tt));
        end
    end
end

% 2) DEFINE THE RELATIVE PERIOD

Tr=ones(size(h0)).*10;

% 3) DEFINE THE ABSOLUTE PERIOD

Ta=Tr;

% 4) COMPUTE THE WAVE NUMBER FOR ZERO CURRENTS

error=1000;
L0=(9.81*Ta.^2)./(2*pi);
L1=L0.*tanh(2*pi*h0./L0);

while(error>0.1)
    L=L0(1,1)*tanh(2*pi*h0(1,1)./L1(1,1));
    error=abs(L-L1(1,1))/L*100;
    L1(1,1)=L;
end

k=zeros(size(u));
k(1,:,:)=(2*pi/L);
k(2,:,:)=(2*pi/L);

wr=zeros(size(u));
wa=zeros(size(u));
wr(1,:,:)=(2*pi)./Tr;
wr(1,:,:)=(2*pi)./Tr(1,1);
wa(1,:,:)=(2*pi)./Tr(1,1);
wr(2,:,:)=(2*pi)./Tr(1,1);

wa=zeros(size(u));
wa=wr;

clear L error

for tt=2:nt-1
    display(tt)
    
    % 5) COMPUTE THE WAVE NUMBER USING NEWTON RAPHSON

    for jj=1:nm
        for ii=1:nn
            error=1000;
            cff_k=k(tt-1,jj,ii);
            cff_wr=wr(tt,jj,ii);
            cff_h=h(tt,jj,ii);
           
            while(error>0.01)
                F=cff_wr.^2-9.81*cff_k*tanh(cff_k*cff_h);
                FD=-9.81*tanh(cff_k*cff_h)-9.81*cff_k*cff_h/(cosh(cff_k*cff_h)).^2;
                cff_k1=cff_k-F/FD;
                error=abs(cff_k1-cff_k)/cff_k*100;
                cff_k=cff_k1;
            end

            k(tt,jj,ii)=cff_k;
        end
    end
    
    % 6) COMPUTE THE WAVE CELERITY
    
    for jj=1:nm
      for ii=1:nn
       c(tt,jj,ii)=9.81*k(tt,jj,ii)*tanh(k(tt,jj,ii)*h(tt,jj,ii))/sinh(2*k(tt,jj,ii)*h(tt,jj,ii));
      end
    end
    

    % 8) COMPUTE THE ABSOLUTE PERIOD TIME VARIATION

    for jj=1:nm
      for ii=1:nn
        wa(tt+1,jj,ii)=wa(tt,jj,ii)+k(tt,jj,ii)*c(tt,jj,ii)*(h(tt+1,jj,ii)-h(tt,jj,ii))...
        +k(tt,jj,ii)*(u(tt+1,jj,ii)-u(tt,jj,ii));
      end
    end    
    
    % 7) COMPUTE THE RELATIVE FREQUENCY

    for jj=1:nm
      for ii=1:nn
        wr(tt,jj,ii)=wa(tt,jj,ii)-k(tt,jj,ii)*u(tt,jj,ii);
      end
    end
    

    % 9) COMPUTE THE RELATIVE PERIOD TIME VARIATION 

    for jj=2:nm
        for ii=2:nn
            wr(tt+1,jj,ii)=wr(tt,jj,ii)-wr(tt,jj,ii)*k(tt,jj,ii).*h(ii,jj,ii)/sinh(2*k(tt,jj,ii).*h(ii,jj,ii))*...
            (((u(tt,jj-1,ii)-u(tt,jj,ii)))*dt/dy+((u(tt,jj,ii-1)-u(tt,jj,ii)))*dt/dx);
        end
    end
    
    wr(tt+1,:,1)=wr(tt+1,:,2);  
    wr(tt+1,1,:)=wr(tt+1,2,:);
    
%     figure
%     pcolor(squeeze(wr(tt+1,:,:)))
%     shading flat
%     title('wr')
%     pause
    
    figure(1)
    
    subplot(2,2,1)
    pcolor(squeeze(h(tt,:,:)))
    shading flat
    colorbar 'vertical'
    caxis([8 12])

    subplot(2,2,2)
    pcolor(squeeze(u(tt,:,:)))
    shading flat
    colorbar 'vertical'
    caxis([-0.2 0.2])
    
    subplot(2,2,3)
    pcolor((2*pi)./squeeze(k(tt,:,:)))
    shading flat
    colorbar 'vertical'
    
    subplot(2,2,4)
    pcolor((2*pi)./squeeze(wr(tt,:,:)))
    shading flat
    colorbar 'vertical'
    caxis([9.5 10.5])
    pause
  
end


figure
plot(2*pi./squeeze(wr(:,31,26)))
hold on
plot(2*pi./squeeze(wa(:,31,26)),'r')
    
