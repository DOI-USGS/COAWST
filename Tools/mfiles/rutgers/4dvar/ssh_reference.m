function [ssh_ref,ssh_err]=ssh_reference(K,rho,ssh,Niter);

%
% SSH_REFERENCE:  Computes the balance operartor reference sea surface height
%
% [ssh_ref,ssh_err]=ssh_reference(K,rho,ssh,Niter)
%
% This function computes the reference sea surface height (m) used in the
% balance operator. It uses ROMS pressure gradient "prsgrd31.h" formulation
% to compute the dynamic pressure from the basic state "in situ" density
% computed by a call to "eos.h" in "ini_balance.h"
%
% On Input:
%
%    K           Balance operator structure array:
%
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
%    rho         In situ density (km/m3), 3D array
%    ssh         First guess sea surface height (m), 2D array
%    Niter       Number of biconjugate gradient iterations, OPTIONAL
%                  (default: K.Niter)
%
% On Output:
%
%    ssh_ref     Reference sea surface height (m), 2D array
%
% Calls:
%
%    biconj:     Biconjugate gradient algorithm
%

% svn $Id: ssh_reference.m 996 2020-01-10 04:28:56Z arango $
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
% its gradient (gradPx) in XI-direction (at U-points).
%
%   (i-1,j,N)      1:L ,1:M,N         DO j=Jstr-1,Jend
%   (i  ,j,N)      2:Lp,1:M,N           DO i=Istr,Iend+1 

cff1(1:L,1:M)=K.Zw(2:Lp,1:M,Np)-K.Zr(2:Lp,1:M,N )+                       ...
              K.Zw(1:L ,1:M,Np)-K.Zr(1:L ,1:M,N );

phix(1:L,1:M)=fac1.*(rho(2:Lp,1:M,N)-                                    ...
                     rho(1:L ,1:M,N)).*cff1(1:L,1:M);

phix_bar(1:L,1:M)=0.5.*(K.Hz(1:L ,1:M,N)+K.Hz(2:Lp,1:M,N)).*             ...
                  phix(1:L,1:M).*                                        ...
                  K.umask(1:L,1:M);

% Compute interior baroclinic pressure gradient (phix, m2/s2) and
% its vertically integrated gradient (gradPx) in the XI-direction
% (at U-points). Differentiate and then vertically integrate.
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

  cff1(1:L,1:M)=(1.0+gamma(1:L,1:M)).*(rho(2:Lp,1:M,k+1)-                ...
                                       rho(1:L ,1:M,k+1))+               ...
                (1.0-gamma(1:L,1:M)).*(rho(2:Lp,1:M,k  )-                ...
                                       rho(1:L ,1:M,k  ));

  cff2(1:L,1:M)=rho(2:Lp,1:M,k+1)+rho(1:L ,1:M,k+1)-                     ...
                rho(2:Lp,1:M,k  )-rho(1:L ,1:M,k  );

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

end,

% Compute surface baroclinic pressure gradient (phie, m2/s2) and
% its gradient (gradPy) in ETA-direction (at V-points).
%
%   (i,j-1,N)      1:L,1:M ,N         DO j=Jstr,Jend+1
%   (i,j  ,N)      1:L,2:Mp,N           DO i=Istr-1,Iend

cff1(1:L,1:M)=K.Zw(1:L,2:Mp,Np)-                                         ...
              K.Zr(1:L,2:Mp,N )+                                         ...
              K.Zw(1:L,1:M ,Np)-                                         ...
              K.Zr(1:L,1:M ,N );
  
phie(1:L,1:M)=fac1.*(rho(1:L,2:Mp,N)-                                    ...
                     rho(1:L,1:M ,N)).*cff1(1:L,1:M);

phie_bar(1:L,1:M)=0.5.*(K.Hz(1:L,1:M ,N)+K.Hz(1:L,2:Mp,N)).*             ...
                  phie(1:L,1:M).*                                        ...
                  K.vmask(1:L,1:M);

% Compute interior baroclinic pressure gradient (phie, m2/s2) and
% its vertically integrated gradient (gradPy) in the ETA-direction
% (at V-points). Differentiate and then vertically integrate.
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

  cff1(1:L,1:M)=(1.0+gamma(1:L,1:M)).*(rho(1:L,2:Mp,k+1)-                ...
                                       rho(1:L,1:M ,k+1))+               ...
                (1.0-gamma(1:L,1:M)).*(rho(1:L,2:Mp,k  )-                ...
                                       rho(1:L,1:M ,k  ));

  cff2(1:L,1:M)=rho(1:L,2:Mp,k+1)+rho(1:L,1:M ,k+1)-                     ...
                rho(1:L,2:Mp,k  )-rho(1:L,1:M ,k  );

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

end,

clear cff1 cff2 cff3 cff4 gamma phie phix

%---------------------------------------------------------------------------
% Compute RHS term (m/s2) for balance sea surface height elliptic equation.
%---------------------------------------------------------------------------
%
% Apply zero boundary conditions.

GradPx(1:L,1:Mp)=0.0;                            % at U-points

GradPx(2:Lm,2:M)=phix_bar(2:Lm,2:M);

GradPy(1:Lp,1:M)=0.0;                            % at V-points

GradPy(2:L,2:Mm)=phie_bar(2:L,2:Mm);

% Computer RHS term (m/s2), at RHO-points.

ssh_rhs(1:Lp,1:Mp)=0.0;

ssh_rhs(2:L,2:M)=-K.pm(2:L,2:M).*K.pn(2:L,2:M).*                         ...
                 (K.pmon_u(2:L ,1:Mm).*GradPx(2:L ,1:Mm)-                ...
                  K.pmon_u(1:Lm,1:Mm).*GradPx(1:Lm,1:Mm)+                ...
                  K.pnom_v(1:Lm,2:M ).*GradPy(1:Lm,2:M )-                ...
                  K.pnom_v(1:Lm,1:Mm).*GradPy(1:Lm,1:Mm)).*              ...
                 K.rmask(2:L,2:M);

%---------------------------------------------------------------------------
% Compute refence sea surface height (m), use biconjugate gradient
% algorithm to solve elliptic equation.
%---------------------------------------------------------------------------

[ssh_ref,ssh_err]=biconj(K,ssh_rhs,ssh,Niter);

return
