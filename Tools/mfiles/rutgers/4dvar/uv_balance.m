function [deltaU_b,deltaV_b,zeta_rhs]=uv_balance(K,deltaR_b)

%
% UV_BALANCE:  Computes the balanced, baroclinic U- and V-momentum anomalies
%
% [deltaU_b,deltaV_b,zeta_rhs]=uv_balance(K,deltaR_b)
%
% This function computes the balanced, baroclinic U- and V-momentum 
% anomalies using geostrophic balance for the error covariance balance
% operator, K^(-1). It uses ROMS pressure gradient "prsgrd31.h" formulation to
% compute the dynamic pressure from the balance density anomaly computed
% in function "rho_balance".
%
% On Input:
%
%    K           Balance operator structure array:
%
%                  K.f         Coriolis parameter (1/s) at RHO-points
%                  K.pm        Curvilinear X-coordinate metric (1/m)
%                  K.pn        Curvilinear Y-coordinate metric (1/m)
%                  K.rmask     Land/Sea mask on RHO-points
%                  K.umask     Land/Sea mask on U-points
%                  K.vmask     Land/Sea mask on V-points
%                  K.Hz        Vertical level thicknesses (m)
%                  K.Zr        Depths at vertical RHO-points (m, negative)
%                  K.Zw        Depths at vertical W-points   (m, negative)
%                  K.g         Accelerarion due to gravity (m/s2)
%                  K.rho0      Mean density (Kg/m3) used when the Boussinesq
%                                approximation is inferred (kg/m3)
%
%    deltaR_b    Balanced density anomaly (kg/m3), as computed from function
%                  "rho_balance", 3D array
%
% On Output:
%
%    deltaU_b    Balanced, baroclinic U-momentum anomaly (m/s), 3D array
%    deltaV_b    Balanced, baroclinic V-momentum anomaly (m/s), 3D array
%    zeta_rhs    RHS term (m4/s2) for balance free-surface elliptic
%                  equation, 2D array
%

% svn $Id: uv_balance.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Compute number of grid points.

[Lp,Mp,Np]=size(K.Zw);

L=Lp-1; Lm=L-1;
M=Mp-1; Mm=M-1;
N=Np-1; Nm=N-1;

% Initialize.

deltaU_b=[];
deltaV_b=[];
zeta_rhs=[];

%---------------------------------------------------------------------------
% Compute pressure gradient components.
%---------------------------------------------------------------------------
%
% NOTE: fac2=0 because the balanced component should consist of the
% baroclinic pressure gradient only.

fac1=0.5*K.g/K.rho0;
fac2=0.0;                        % originally, fac2=g (not used anyway)
fac3=0.25*K.g/K.rho0;

% Compute surface baroclinic pressure gradient (phix, m2/s2) and
% its vertically integrated values (phix_bar, m3/s2; gradPx, m/s)
% in the XI-direction (at U-points).
%
%   (i-1,j,N)      1:L ,1:M,N         DO j=Jstr-1,Jend
%   (i  ,j,N)      2:Lp,1:M,N           DO i=Istr,Iend+1 

cff1(1:L,1:M)=K.Zw(2:Lp,1:M,Np)-                                         ...
              K.Zr(2:Lp,1:M,N )+                                         ...
              K.Zw(1:L ,1:M,Np)-                                         ...
              K.Zr(1:L ,1:M,N );

phix(1:L,1:M)=fac1.*(deltaR_b(2:Lp,1:M,N)-                               ...
                     deltaR_b(1:L ,1:M,N)).*cff1(1:L,1:M);

phix_bar(1:L,1:M)=0.5.*(K.Hz(1:L ,1:M,N)+K.Hz(2:Lp,1:M,N)).*             ...
                  phix(1:L,1:M).*                                        ...
                  K.umask(1:L,1:M);

gradPx(1:L,1:M,N)=0.5.*phix(1:L,1:M).*                                   ...
                       (K.pm(1:L ,1:M)+K.pm(2:Lp,1:M))./                 ...
                       ( K.f(1:L ,1:M)+ K.f(2:Lp,1:M));

% Compute interior baroclinic pressure gradient (phix, m2/s2) and
% its vertically integrated values (phix_bar, m3/s2; gradPx, m/s)
% in the XI-direction (at U-points). Differentiate and then
% vertically integrate.
%
%   (i-1,j,k  )    1:L ,1:M,k        DO k=1,N-1
%   (i  ,j,k  )    2:Lp,1:M,k          DO j=Jstr-1,Jend
%   (i-1,j,k+1)    1:L ,1:M,k+1          DO i=Istr,Iend+1
%   (i  ,j,k+1)    2:Lp,1:M,k+1

for k=Nm:-1:1,

  cff1(1:L,1:M)=1.0./((K.Zr(2:Lp,1:M,k+1)-K.Zr(2:Lp,1:M,k  )).*          ...
                      (K.Zr(1:L ,1:M,k+1)-K.Zr(1:L ,1:M,k  )));

  cff2(1:L,1:M)=K.Zr(2:Lp,1:M,k  )-K.Zr(1:L ,1:M,k  )+                   ...
                K.Zr(2:Lp,1:M,k+1)-K.Zr(1:L ,1:M,k+1);

  cff3(1:L,1:M)=K.Zr(2:Lp,1:M,k+1)-K.Zr(2:Lp,1:M,k  )-                   ...
                K.Zr(1:L ,1:M,k+1)+K.Zr(1:L ,1:M,k  );

  gamma(1:L,1:M)=0.125.*cff1.*cff2.*cff3;

  cff1(1:L,1:M)=(1.0+gamma(1:L,1:M)).*(deltaR_b(2:Lp,1:M,k+1)-           ...
                                       deltaR_b(1:L ,1:M,k+1))+          ...
                (1.0-gamma(1:L,1:M)).*(deltaR_b(2:Lp,1:M,k  )-           ...
                                       deltaR_b(1:L ,1:M,k  ));

  cff2(1:L,1:M)=deltaR_b(2:Lp,1:M,k+1)+deltaR_b(1:L ,1:M,k+1)-           ...
                deltaR_b(2:Lp,1:M,k  )-deltaR_b(1:L ,1:M,k  );

  cff3(1:L,1:M)=K.Zr(2:Lp,1:M,k+1)+K.Zr(1:L ,1:M,k+1)-                   ...
                K.Zr(2:Lp,1:M,k  )-K.Zr(1:L ,1:M,k  );

  cff4(1:L,1:M)=(1.0+gamma(1:L,1:M)).*(K.Zr(2:Lp,1:M,k+1)-               ...
                                       K.Zr(1:L ,1:M,k+1))+              ...
                (1.0-gamma(1:L,1:M)).*(K.Zr(2:Lp,1:M,k  )-               ...
                                       K.Zr(1:L ,1:M,k  ));

  phix(1:L,1:M)=phix(1:L,1:M)+fac3.*(cff1.*cff3-cff2.*cff4);

  phix_bar(1:L,1:M)=phix_bar(1:L,1:M)+                                   ...
                    0.5.*(K.Hz(1:L ,1:M,k)+K.Hz(2:Lp,1:M,k)).*           ...
                    phix(1:L,1:M).*                                      ...
                    K.umask(1:L,1:M);

  gradPx(1:L,1:M,k)=0.5.*phix(1:L,1:M).*                                 ...
                         (K.pm(1:L ,1:M)+K.pm(2:Lp,1:M))./               ...
                         ( K.f(1:L ,1:M)+ K.f(2:Lp,1:M));

end,

% Compute surface baroclinic pressure gradient (phie, m2/s2) and
% its vertically integrated value (phie_bar, m3/s2; gradPy, m/s)
% in the ETA-direction (at V-points).
%
%   (i,j-1,N)      1:L,1:M ,N         DO j=Jstr,Jend+1
%   (i,j  ,N)      1:L,2:Mp,N           DO i=Istr-1,Iend

cff1(1:L,1:M)=K.Zw(1:L,2:Mp,Np)-                                         ...
              K.Zr(1:L,2:Mp,N )+                                         ...
              K.Zw(1:L,1:M ,Np)-                                         ...
              K.Zr(1:L,1:M ,N );
  
phie(1:L,1:M)=fac1.*(deltaR_b(1:L,2:Mp,N)-                               ...
                     deltaR_b(1:L,1:M ,N)).*cff1(1:L,1:M);

phie_bar(1:L,1:M)=0.5.*(K.Hz(1:L,1:M ,N)+K.Hz(1:L,2:Mp,N)).*             ...
                  phie(1:L,1:M).*                                        ...
                  K.vmask(1:L,1:M);

gradPy(1:L,1:M,N)=0.5.*phie(1:L,1:M).*                                   ...
                       (K.pn(1:L,1:M )+K.pn(1:L,2:Mp))./                 ...
                       ( K.f(1:L,1:M )+ K.f(1:L,2:Mp));

% Compute interior baroclinic pressure gradient (phie, m2/s2) and
% its vertically integrated value (phie_bar, m3/s2; gradPy, m/s)
% in the ETA-direction (at V-points). Differentiate and then
% vertically integrate.
%
%   (i,j-1,k  )    1:L,1:M ,k         DO k=1,N-1
%   (i,j  ,k  )    1:L,2:Mp,k           DO j=Jstr,Jend+1
%   (i,j-1,k+1)    1:L,1:M ,k+1           DO i=Istr-1,Iend
%   (i,j  ,k+1)    1:L,2:Mp,k+1

for k=Nm:-1:1,

  cff1(1:L,1:M)=1.0./((K.Zr(1:L,2:Mp,k+1)-K.Zr(1:L,2:Mp,k  )).*          ...
                      (K.Zr(1:L,1:M ,k+1)-K.Zr(1:L,1:M ,k  )));
    
  cff2(1:L,1:M)=K.Zr(1:L,2:Mp,k  )-K.Zr(1:L,1:M ,k  )+                   ...
                K.Zr(1:L,2:Mp,k+1)-K.Zr(1:L,1:M ,k+1);

  cff3(1:L,1:M)=K.Zr(1:L,2:Mp,k+1)-K.Zr(1:L,2:Mp,k  )-                   ...
                K.Zr(1:L,1:M ,k+1)+K.Zr(1:L,1:M ,k  );
   
  gamma(1:L,1:M)=0.125.*cff1.*cff2.*cff3;

  cff1(1:L,1:M)=(1.0+gamma(1:L,1:M)).*(deltaR_b(1:L,2:Mp,k+1)-           ...
                                       deltaR_b(1:L,1:M ,k+1))+          ...
                (1.0-gamma(1:L,1:M)).*(deltaR_b(1:L,2:Mp,k  )-           ...
                                       deltaR_b(1:L,1:M ,k  ));

  cff2(1:L,1:M)=deltaR_b(1:L,2:Mp,k+1)+deltaR_b(1:L,1:M ,k+1)-           ...
                deltaR_b(1:L,2:Mp,k  )-deltaR_b(1:L,1:M ,k  );

  cff3(1:L,1:M)=K.Zr(1:L,2:Mp,k+1)+K.Zr(1:L,1:M ,k+1)-                   ...
                K.Zr(1:L,2:Mp,k  )-K.Zr(1:L,1:M ,k  );

  cff4(1:L,1:M)=(1.0+gamma(1:L,1:M)).*(K.Zr(1:L,2:Mp,k+1)-               ...
                                       K.Zr(1:L,1:M ,k+1))+              ...
                (1.0-gamma(1:L,1:M)).*(K.Zr(1:L,2:Mp,k  )-               ...
                                       K.Zr(1:L,1:M ,k  ));

  phie(1:L,1:M)=phie(1:L,1:M)+                                           ...
                fac3.*(cff1.*cff3-cff2.*cff4);

  phie_bar(1:L,1:M)=phie_bar(1:L,1:M)+                                   ...
                    0.5.*(K.Hz(1:L,1:M ,k)+K.Hz(1:L,2:Mp,k)).*           ...
                    phie(1:L,1:M).*                                      ...
                    K.vmask(1:L,1:M);

  gradPy(1:L,1:M,k)=0.5.*phie(1:L,1:M).*                                 ...
                         (K.pn(1:L,1:M )+K.pn(1:L,2:Mp))./               ...
                         ( K.f(1:L,1:M )+ K.f(1:L,2:Mp));

end,

clear cff1 cff2 cff3 cff4 gamma phie phix

%---------------------------------------------------------------------------
% Compute balance horizontal momentum anomalies.
%---------------------------------------------------------------------------
%
% Set zero boundary conditions.

deltaU_b(1:L,1:Mp,1:N)=0.0;

mask=repmat(K.umask,[1,1,N]);

deltaU_b(2:Lm,2:M,1:N)=-0.25.*(gradPy(1:Lm-1,1:Mm,1:N)+                  ...
                               gradPy(2:Lm  ,1:Mm,1:N)+                  ...
                               gradPy(1:Lm-1,2:M ,1:N)+                  ...
                               gradPy(2:Lm  ,2:M ,1:N)).*                ...
                              mask(2:Lm,2:M,1:N);

deltaV_b(1:Lp,1:M,1:N)=0.0;

mask=repmat(K.vmask,[1,1,N]);

deltaV_b(2:L,2:Mm,1:N)=0.25.*(gradPx(1:Lm,1:Mm-1,1:N)+                   ...
                              gradPx(2:L ,1:Mm-1,1:N)+                   ...
                              gradPx(1:Lm,2:Mm  ,1:N)+                   ...
                              gradPx(2:L ,2:Mm  ,1:N)).*                 ...
                             mask(2:L,2:Mm,1:N);

clear mask gradPx gradPy

%---------------------------------------------------------------------------
% Compute RHS term (m/s2) for balance sea surface height elliptic equation:
%
%     div (h grad(ssh)) = - div (phi_bar)
%
% where phi_bar is the vertically integrated dynamic pressure vector.
%---------------------------------------------------------------------------
%
% Apply zero boundary conditions.

GradPx(1:L,1:Mp)=0.0;                            % at U-points

GradPx(2:Lm,2:M)=phix_bar(2:Lm,2:M);

GradPy(1:Lp,1:M)=0.0;                            % at V-points

GradPy(2:L,2:Mm)=phie_bar(2:L,2:Mm);

% Computer RHS term (m/s2), at RHO-points.

zeta_rhs(1:Lp,1:Mp)=0.0;

zeta_rhs(2:L,2:M)=-K.pm(2:L,2:M).*K.pn(2:L,2:M).*                        ...
                  (K.pmon_u(2:L ,1:Mm).*GradPx(2:L ,1:Mm)-               ...
                   K.pmon_u(1:Lm,1:Mm).*GradPx(1:Lm,1:Mm)+               ...
                   K.pnom_v(1:Lm,2:M ).*GradPy(1:Lm,2:M )-               ...
                   K.pnom_v(1:Lm,1:Mm).*GradPy(1:Lm,1:Mm)).*             ...
                  K.rmask(2:L,2:M);

return
