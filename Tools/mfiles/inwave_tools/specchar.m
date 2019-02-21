function [Hs,Tp]=specchar(E2,E1,df,dtheta,nfreq,ndir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    CALCULATE MEAN WAVE PARAMETERS   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Compute Hs and Mwd from 2-dimensional spectrum
Etotal=0;
Mwd_C=0;
Mwd_S=0;
for j=1:ndir
    for i=1:nfreq
        Etotal=Etotal+(E2(i,j))*dtheta(j)*df(i);
    end
end
%Hs=4*sqrt(Etotal/1025/9.81)
Hs=4*sqrt(Etotal)

%%% Compute Tp from 1-dimensional spectrum
[Emax Emax_pos]=max(E1);
Tp=1./E1(Emax_pos(2),1)
end