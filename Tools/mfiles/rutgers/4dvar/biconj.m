function [ssh,err]=biconj(K,rhs_r2d,ssh_guess,Niter);

%
% BICONJ:  Biconjugate gradient solver for the SSH elliptic equation
% 
% This function uses a bicojugate gradient algorithm to solve the 
% sea surface height elliptic equation (Fukumori et al., 1998):
%
%      div (h grad(ssh)) = - div (int{int{grad(rho)z'/rho0) dz'} dz})
%
% On Input:
%
%    K         Balance operator structure array:
%
%              K.Lm           Number of interior points in XI-direction
%              K.Mm           Number of interior points in ETA-direction
%              K.EW_periodic  Switch for East-West periodic boundaries
%              K.NS_periodic  Switch for North-South periodic boundaries
%              K.g            Accelerarion due to gravity (m/s2)
%              K.h            Bathymetry (m) at RHO-points
%              K.pm           Curvilinear X-coordinate metric (1/m)
%              K.pn           Curvilinear Y-coordinate metric (1/m)
%              K.pmon_u       Curvilinear metric pm/pn at U-points
%              K.pnom_v       Curvilinear metric pn/pm at V-points
%              K.rmask        Land/Sea mask on RHO-points
%              K.umask        Land/Sea mask on U-points
%              K.vmask        Land/Sea mask on V-points
%
%    rhs_r2d   RHS term for elliptic equation (m3/s2) as computed
%                from the pressure gradient in function "uv_balance":
%
%                int{int{grad(rho)z'/rho0) dz'} dz}
%
%    ssh_guess Sea surface height (m) first guess, use backgound as
%                starting value (2D array)
%    Niter     Number of biconjugate gradient iterations, OPTIONAL
%                (default: K.Niter)
%
% On Output:
%
%    ssh       Baroclinic sea surface height (m), 2D array
%    err       SSH error after Niter iterations (scalar)
%
% Reference:
%
%   Fukumori, I., R. Raghunath and L. Fu, 1998: Nature of global
%     large-scale sea level variability in relation to atmospheric
%     forcing: a modeling study, J. Geophys. Res., 103, 5493-5512.
%
% Internal routines elliptical solver:
%
%    r2d_oper:   Computes SSH elliptic equation operator and
%                   its transpose.
%
%    r2d_bc:     Sets boundary conditions
%
%    r2d_dotp    Computes dot products for biconjugate gradient
%                  algorithm.
%
%    ad_r2d_bc:  Sets adjoint boundary conditions
%

% svn $Id: biconj.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license           Andrew M. Moore         %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Initialize internal parameters.

if (nargin < 4),
  Niter=K.Niter;
end,

Lm=K.Lm; L=Lm+1; Lp=L+1;
Mm=K.Mm; M=Mm+1; Mp=M+1;

%  Initialize local arrays.

v_r2d=zeros(size(ssh_guess));
bv_r2d=zeros(size(ssh_guess));
z1_r2d=zeros(size(ssh_guess));
z2_r2d=zeros(size(ssh_guess));
z3_r2d=zeros(size(ssh_guess));

%  Choose the starting value (m) of "r2d_ref". Use first guess sea surface
%  height.

r2d_ref(2:L,2:M)=ssh_guess(2:L,2:M).*K.rmask(2:L,2:M);

[r2d_ref]=r2d_bc(r2d_ref,2,L,2,M,K.EW_periodic,K.NS_periodic);

%  Compute starting value (m/s2) of divergence operator: 
%     z1_r2d = div[h grad(r2d_ref)].

Ltrans=0;
[pc_r2d,z1_r2d]=r2d_oper(K,Ltrans,r2d_ref);

%  Set the initial values for residual vectors "r" and "br". Then, use
%  recurrence relationship to compute direction vectors "p" and "bp".

r_r2d (2:L,2:M,1)=rhs_r2d(2:L,2:M)-z1_r2d(2:L,2:M);

br_r2d(2:L,2:M,1)=r_r2d (2:L,2:M,1);

p_r2d (2:L,2:M,1)=r_r2d (2:L,2:M,1)./pc_r2d(2:L,2:M);

bp_r2d(2:L,2:M,1)=br_r2d(2:L,2:M,1)./pc_r2d(2:L,2:M);

%===========================================================================
%  Iterate.
%===========================================================================

for iter=1:Niter-1,

  z1_r2d(2:L,2:M)=p_r2d(2:L,2:M,iter);

  [z1_r2d]=r2d_bc(z1_r2d,2,L,2,M,K.EW_periodic,K.NS_periodic);

%  Compute divergence operator: v_r2d = div[h grad(z1_r2d)].

  Ltrans=0;
  [pc_r2d,v_r2d]=r2d_oper(K,Ltrans,z1_r2d);
   
%  Compute dot products and "bc_ak" coefficient.

  z1_r2d(2:L,2:M)=r_r2d (2:L,2:M,iter)./pc_r2d(2:L,2:M);

  z2_r2d(2:L,2:M)=br_r2d(2:L,2:M,iter);
   
  z3_r2d(2:L,2:M)=bp_r2d(2:L,2:M,iter);

  [zdf1(iter)]=r2d_dotp(K,z2_r2d,z1_r2d);

  [zdf2(iter)]=r2d_dotp(K,z3_r2d, v_r2d);

  bc_ak(iter)=zdf1(iter)/zdf2(iter);

%  Solve for new iterate of "r2d_ref".

  r2d_ref(2:L,2:M)=r2d_ref(2:L,2:M)+                                     ...
                   bc_ak(iter).*p_r2d(2:L,2:M,iter);

%---------------------------------------------------------------------------
%  Test for convergence: use "bv_r2d" as temporary storage.
%---------------------------------------------------------------------------

  if (iter == Niter-1),

    [r2d_ref]=r2d_bc(r2d_ref,2,L,2,M,K.EW_periodic,K.NS_periodic);
      
%  Compute divergence operator: bv_r2d = div[h grad(r2d_ref)].

    Ltrans=0;
    [pc_r2d,bv_r2d]=r2d_oper(K,Ltrans,r2d_ref);

%  Compute dot products and report convergence value.

    bv_r2d(2:L,2:M)=bv_r2d(2:L,2:M)-rhs_r2d(2:L,2:M);

    [zdf4]=r2d_dotp(K,bv_r2d,bv_r2d);

    [zdf5]=r2d_dotp(K,rhs_r2d,rhs_r2d);

    err=sqrt(zdf4/zdf5);
       
    disp(['   BICONJ - Error in reference sea surface height = ',        ...
          sprintf('%15.8e', err),                                        ...
          ',  Rec = ', sprintf('%3.3i', K.TimeRec), ]);

  end,

%---------------------------------------------------------------------------
%  Compute new (iter+1) residual and direction vectors.
%---------------------------------------------------------------------------

  z1_r2d(2:L,2:M)=bp_r2d(2:L,2:M,iter);

  [z1_r2d]=r2d_bc(z1_r2d,2,L,2,M,K.EW_periodic,K.NS_periodic);

%  Compute divergence operator, tranpose: bv_r2d = div[h grad(z1_r2d)].
%  Need to call "ad_r2d_bc" here since Ltrans is TRUE.

  Ltrans=1;
  [pc_r2d,bv_r2d]=r2d_oper(K,Ltrans,z1_r2d);
     
%  Compute new residual vectors "r" and "br".

  r_r2d (2:L,2:M,iter+1)=r_r2d (2:L,2:M,iter)-                           ...
                         bc_ak(iter).*v_r2d (2:L,2:M);

  br_r2d(2:L,2:M,iter+1)=br_r2d(2:L,2:M,iter)-                           ...
                         bc_ak(iter).*bv_r2d(2:L,2:M);

%  Compute dot product and "bc_ak" coefficient.

  z1_r2d(2:L,2:M)=r_r2d (2:L,2:M,iter+1)./pc_r2d(2:L,2:M);

  z2_r2d(2:L,2:M)=br_r2d(2:L,2:M,iter+1);

  [zdf3(iter)]=r2d_dotp(K,z1_r2d,z2_r2d);
     
  bc_bk(iter)=zdf3(iter)./zdf1(iter);

%  Use recurrence relationship to compute new direction vectors
%  "p" and "bp". 

  p_r2d (2:L,2:M,iter+1)=r_r2d (2:L,2:M,iter+1)./pc_r2d(2:L,2:M)+        ...
                         bc_bk(iter).*p_r2d (2:L,2:M,iter);

  bp_r2d(2:L,2:M,iter+1)=br_r2d(2:L,2:M,iter+1)./pc_r2d(2:L,2:M)+        ...
                         bc_bk(iter).*bp_r2d(2:L,2:M,iter);

end, 

%---------------------------------------------------------------------------
%  Set baroclinic sea surface height (m);
%---------------------------------------------------------------------------

ssh=r2d_ref;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pc_r2d,r2d_out]=r2d_oper(K,Ltrans,r2d_in);

%
% This function computes SSH elliptic equation divergence operator
%
%              div(h grad(r2d_in))
%
% and its transpose:
%
% On Input:
%
%    K         Balance operator structure array:
%
%              K.Lm           Number of interior points in the XI-direction
%              K.Mm           Number of interior points in the ETA-direction
%              K.g            Accelerarion due to gravity (m/s2)
%              K.pm           Curvilinear X-coordinate metric (1/m)
%              K.pn           Curvilinear Y-coordinate metric (1/m)
%              K.pmon_u       Curvilinear metric pm/pn at U-points
%              K.pnom_v       Curvilinear metric pn/pm at V-points
%              K.rmask        Land/Sea mask on RHO-points
%              K.umask        Land/Sea mask on U-points
%              K.vmask        Land/Sea mask on V-points
%
%    Ltrans    Switch for transposed divergence operator (modified adjoint)
%    r2d_in    Current SSH iterate (m), 2D array
%
% On Output:
%
%    pc_r2d    Elliptic equation scale coefficient (1/s2), 2D array
%    r2d_out   SSH divergence (m/s2), 2D array
%

%  Initialize.

cff=0.5*K.g;

Lm=K.Lm; L=Lm+1; Lp=L+1;
Mm=K.Mm; M=Mm+1; Mp=M+1;

r2d_out=zeros(size(r2d_in));
FX=zeros(size(r2d_in));
FE=zeros(size(r2d_in));

%---------------------------------------------------------------------------
%  Compute divergence XI- and ETA-components: transposed operator.
%                                             (Ltrans=1)
%---------------------------------------------------------------------------

if (Ltrans),

%> for j=Jstr:Jend,
%>   for i=Istr:Iend,
%>     r2d_out(i,j)=K.pm(i,j)*K.pn(i,j)*                                 ...
%>                  (FX(i+1,j)-FX(i,j)+                                  ...
%>                   FE(i,j+1)-FE(i,j))*                                 ...
%>                  K.rmask(i,j);
%>   end,
%> end,
%>                                           adjoint code below

  fac(1:Lm,1:Mm)=K.pm(2:L,2:M).*                                         ...
                 K.pn(2:L,2:M).*                                         ...
                 r2d_in(2:L,2:M).*                                       ...
                 K.rmask(2:L,2:M);

  FX(2:L ,2:M )=FX (2:L ,2:M )-                                          ...
                fac(1:Lm,1:Mm);

  FX(3:Lp,2:M )=FX (3:Lp,2:M )+                                          ...
                fac(1:Lm,1:Mm);

  FE(2:L ,2:M )=FE (2:L ,2:M )-                                          ...
                fac(1:Lm,1:Mm);

  FE(2:L ,3:Mp)=FE (2:L ,3:Mp)+                                          ...
                fac(1:Lm,1:Mm);

%> for j=Jstr:Jend+1,
%>   for i=Istr,Iend,
%>     FE(i,j)=cff*K.pnom_v(i,j)*(K.h(i,j)+K.h(i,j-1))*                  ...
%>             (r2d_in(i,j)-r2d_in(i,j-1))*                              ...
%>             vmask(i,j);
%>   end,
%> end,
%>                                           adjoint code below

  fac(1:Lm,1:M)=cff.*K.pnom_v(2:L,1:M).*                                 ...
                (K.h(2:L,2:Mp)+                                          ...
                 K.h(2:L,1:M )).*                                        ...
                 FE(2:L,2:Mp).*                                          ...
                K.vmask(2:L,1:M);

  r2d_out(2:L,1:M )=r2d_out(2:L ,1:M )-                                  ...
                        fac(1:Lm,1:M );

  r2d_out(2:L,2:Mp)=r2d_out(2:L ,2:Mp)+                                  ...
                        fac(1:Lm,1:M );

%> for j=Jstr:Jend,
%>   for i=Istr:Iend+1,
%>     FX(i,j)=cff*K.pmon_u(i,j)*(K.h(i,j)+K.h(i-1,j))*                  ...
%>             (r2d_in(i,j)-r2d_in(i-1,j))*                              ...
%>             K.umask(i,j);
%>   end,
%> end,
%>                                           adjoint code below

  fac(1:L,1:Mm)=cff.*K.pmon_u(1:L,2:M).*                                 ...
                (K.h(2:Lp,2:M)+                                          ...
                 K.h(1:L ,2:M)).*                                        ...
                FX(2:Lp,2:M).*                                           ...
                K.umask(1:L,2:M);

  r2d_out(2:Lp,2:M)=r2d_out(2:Lp,2:M )+                                  ...
                        fac(1:L ,1:Mm);

  r2d_out(1:L ,2:M)=r2d_out(1:L ,2:M )-                                  ...
                        fac(1:L ,1:Mm);

%---------------------------------------------------------------------------
%  Compute divergence XI- and ETA-components: regular operator
%                                             (Ltrans=0).
%---------------------------------------------------------------------------

else,      

%  for j=Jstr:Jend,
%    for i=Istr:Iend+1,
%      FX(i,j)=cff.*K.pmon_u(i,j).*(K.h(i,j)+K.h(i-1,j)).*               ...
%              (r2d_in(i,j)-r2d_in(i-1,j)).*                             ...
%              K.umask(i,j);
%    end,
%  end,
%                                            Fortran like code for reference

  FX(2:Lp,2:M)=cff.*K.pmon_u(1:L,2:M).*                                  ...
               (K.h(2:Lp,2:M)+                                           ...
                K.h(1:L ,2:M)).*                                         ...
               (r2d_in(2:Lp,2:M)-                                        ...
                r2d_in(1:L ,2:M)).*                                      ...
               K.umask(1:L,2:M);

%  for j=Jstr:Jend+1,
%    for i=Istr:Iend,
%      FE(i,j)=cff.*K.pnom_v(i,j).*(K.h(i,j)+K.h(i,j-1)).*               ...
%              (r2d_in(i,j)-r2d_in(i,j-1)).*                             ...
%              K.vmask(i,j);
%    end,
%  end,
%                                            Fortran like code for reference

  FE(2:L,2:Mp)=cff*K.pnom_v(2:L,1:M).*                                   ...
               (K.h(2:L,2:Mp)+                                           ...
                K.h(2:L,1:M )).*                                         ...
               (r2d_in(2:L,2:Mp)-                                        ...
                r2d_in(2:L,1:M )).*                                      ...
               K.vmask(2:L,1:M);

%  Compute divergence (m/s2).

%  for j=Jstr:Jend,
%    for i=Istr:Iend,
%      r2d_out(i,j)=K.pm(i,j).*K.pn(i,j).*                               ...
%                   (FX(i+1,j)-FX(i,j)+                                  ...
%                    FE(i,j+1)-FE(i,j)).*                                ...
%                   K.rmask(i,j);
%    end,
%  end,
%                                            Fortran like code for reference

  r2d_out(2:L,2:M)=K.pm(2:L,2:M).*                                       ...
                   K.pn(2:L,2:M).*                                       ...
                   (FX(3:Lp,2:M )-                                       ...
                    FX(2:L ,2:M )+                                       ...
                    FE(2:L ,3:Mp)-                                       ...
                    FE(2:L ,2:M)).*                                      ...
                   K.rmask(2:L,2:M);
	     
end,

%---------------------------------------------------------------------------
%  Compute scale (1/s2) coefficient.
%---------------------------------------------------------------------------

%  for j=Jstr:Jend,
%    for i=Istr:Iend,
%      pc_r2d(i,j)=-K.pm(i,j).*K.pn(i,j)*                                ...
%                  cff*(K.pnom_v(i  ,j+1)*(K.h(i  ,j+1)+K.h(i  ,j  ))+   ...
%                       K.pnom_v(i  ,j  )*(K.h(i  ,j  )+K.h(i  ,j-1))+   ...
%                       K.pmon_u(i+1,j  )*(K.h(i+1,j  )+K.h(i  ,j  ))+   ...
%                       K.pmon_u(i  ,j  )*(K.h(i  ,j  )+K.h(i-1,j  )));
%    end,
%  end,
%                                            Fortran like code for reference

pc_r2d(2:L,2:M)=-K.pm(2:L,2:M).*                                         ...
                 K.pn(2:L,2:M).*                                         ...
                cff.*(K.pnom_v(2:L,2:M ).*                               ...
                      (K.h(2:L,3:Mp)+                                    ...
                       K.h(2:L,2:M ))+                                   ...
                      K.pnom_v(2:L,1:Mm).*                               ...
                      (K.h(2:L,2:M )+                                    ...
                       K.h(2:L,1:Mm))+                                   ...
                      K.pmon_u(2:L ,2:M).*                               ...
                      (K.h(3:Lp,2:M)+                                    ...
                       K.h(2:L ,2:M))+                                   ...
                      K.pmon_u(1:Lm,2:M).*                               ...
                      (K.h(2:L ,2:M)+                                    ...
                       K.h(1:Lm,2:M)));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A]=r2d_bc(A,Istr,Iend,Jstr,Jend,EW_periodic,NS_periodic);

%
% This function applies boundary conditions for 2D elliptic operator.
%
% On Input:
%
%    A             Field to process (2D array)
%    Istr          Starting interior point I-index (integer)
%    Iend          Ending   interior point I-index (integer)
%    Jstr          Starting interior point J-index (integer)
%    Jend          Ending   interior point J-index (integer)
%    EW_periodic   Switch to apply E-W periodic boundary conditions
%                    (OPTIONAL)
%    NS_periodic   Switch to apply N-S periodic boundary conditions
%                    (OPTIONAL)
%
% On Output:
%
%    A             Field with updated boundary conditions (2D array)
%

%  If not provided, turn off periodic boundary conditions.

if (nargin < 7),
  NS_periodic=0;
end,

if (nargin < 6),
  EW_periodic=0;
end,

%---------------------------------------------------------------------------
%  East-West boundary conditions.
%---------------------------------------------------------------------------

if (~ EW_periodic),
  A(Iend+1,Jstr:Jend)=0.0;
  A(Istr-1,Jstr:Jend)=0.0;
end,

%---------------------------------------------------------------------------
%  North-South boundary conditions.
%---------------------------------------------------------------------------

if (~ NS_periodic),
  A(Istr:Iend,Jend+1)=0.0;
  A(Istr:Iend,Jstr-1)=0.0;
end,

%---------------------------------------------------------------------------
%  Boundary coorners.
%---------------------------------------------------------------------------

if (~ EW_periodic & ~ NS_periodic),
  A(Istr-1,Jstr-1)=0.5.*(A(Istr  ,Jstr-1)+A(Istr-1,Jstr  ));
  A(Iend+1,Jstr-1)=0.5.*(A(Iend  ,Jstr-1)+A(Iend+1,Jstr  ));
  A(Istr-1,Jend+1)=0.5.*(A(Istr-1,Jend  )+A(Istr  ,Jend+1));
  A(Iend+1,Jend+1)=0.5.*(A(Iend+1,Jend  )+A(Iend  ,Jend+1));
end,

%---------------------------------------------------------------------------
%  Periodic boundary conditions.
%---------------------------------------------------------------------------

if (EW_periodic | NS_periodic),
% [A]=periodic(A,Istr,Iend,Jstr,Jend,EW_periodic,NS_periodic);
end,

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ad_A]=ad_r2d_bc(ad_A,Istr,Iend,Jstr,Jend,EW_periodic,NS_periodic);

%
% This function applies adjoint boundary conditions for 2D elliptic
% operator.
%
% On Input:
%
%    ad_A          Adjoint field to process (2D array)
%    Istr          Starting interior point I-index (integer)
%    Iend          Ending   interior point I-index (integer)
%    Jstr          Starting interior point J-index (integer)
%    Jend          Ending   interior point J-index (integer)
%    EW_periodic   Switch to apply E-W periodic boundary conditions
%                    (OPTIONAL)
%    NS_periodic   Switch to apply N-S periodic boundary conditions
%                    (OPTIONAL)
%
% On Output:
%
%    ad_A          Field with updated adjoint boundary conditions
%                    (2D array)
%

%  If not provided, turn off periodic boundary conditions.

if (nargin < 7),
  NS_periodic=0;
end,

if (nargin < 6),
  EW_periodic=0;
end,

%---------------------------------------------------------------------------
%  Adjoint periodic boundary conditions.
%---------------------------------------------------------------------------

if (EW_periodic | NS_periodic),

%>  [A]=periodic(A,Istr,Iend,Jstr,Jend,EW_periodic,NS_periodic);
%>

% [ad_A]=periodic(ad_A,Istr,Iend,Jstr,Jend,EW_periodic,NS_periodic);

end,

%---------------------------------------------------------------------------
%  Boundary corners.
%---------------------------------------------------------------------------

if (~ EW_periodic & ~ NS_periodic),

%>  A(Iend+1,Jend+1)=0.5.*(A(Iend+1,Jend  )+A(Iend  ,Jend+1));
%>

  adfac=0.5*ad_A(Iend+1,Jend+1);
  ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac;
  ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac;
  ad_A(Iend+1,Jend+1)=0.0;

%>  A(Istr-1,Jend+1)=0.5.*(A(Istr-1,Jend  )+A(Istr  ,Jend+1));
%>

  adfac=0.5*ad_A(Istr-1,Jend+1);
  ad_A(Istr-1,Jend  )=ad_A(Istr-1,Jend  )+adfac;
  ad_A(Istr  ,Jend+1)=ad_A(Istr  ,Jend+1)+adfac;
  ad_A(Istr-1,Jend+1)=0.0;

%>  A(Iend+1,Jstr-1)=0.5.*(A(Iend  ,Jstr-1)+A(Iend+1,Jstr  ));
%>

  adfac=0.5*ad_A(Iend+1,Jstr-1);
  ad_A(Iend  ,Jstr-1)=ad_A(Iend  ,Jstr-1)+adfac;
  ad_A(Iend+1,Jstr  )=ad_A(Iend+1,Jstr  )+adfac;
  ad_A(Iend+1,Jstr-1)=0.0;

%>  A(Istr-1,Jstr-1)=0.5.*(A(Istr  ,Jstr-1)+A(Istr-1,Jstr  ));
%>

  adfac=0.5*ad_A(Istr-1,Jstr-1);
  ad_A(Istr  ,Jstr-1)=ad_A(Istr  ,Jstr-1)+adfac;
  ad_A(Istr-1,Jstr  )=ad_A(Istr-1,Jstr  )+adfac;
  ad_A(Istr-1,Jstr-1)=0.0;

end,

%---------------------------------------------------------------------------
%  North-South boundary conditions.
%---------------------------------------------------------------------------

if (~ NS_periodic),

%>  A(Istr:Iend,Jstr-1)=0.0;
%>

  ad_A(Istr:Iend,Jstr-1)=0.0;
  
%>  A(Istr:Iend,Jend+1)=0.0;
%>

  ad_A(Istr:Iend,Jend+1)=0.0;

end,

%---------------------------------------------------------------------------
%  East-West boundary conditions.
%---------------------------------------------------------------------------

if (~ EW_periodic),

%>  A(Istr-1,Jstr:Jend)=0.0;
%>

  ad_A(Istr-1,Jstr:Jend)=0.0;

%>  A(Iend+1,Jstr:Jend)=0.0;
%>

  ad_A(Iend+1,Jstr:Jend)=0.0;

end,

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DotProd]=r2d_dotp(K,s1,s2);

%
% This function computes the dot product betwen s1 and s2.
%
% On Input:
%
%    K         Balance operator structure array:
%
%              K.Lm           Number of interior points in XI-direction
%              K.Mm           Number of interior points in ETA-direction
%              K.rmask        Land/Sea mask on RHO-points
%
%    s1        Biconjugate gradient direction (2D array)
%    s2        Biconjugate gradient direction (2D array)
%
% On Output:
%
%    DotProd   Dot product between s1 and s2.
%

%  Initialize internal parameters.

Lm=K.Lm; L=Lm+1; Lp=L+1;
Mm=K.Mm; M=Mm+1; Mp=M+1;

%---------------------------------------------------------------------------
%  Compute dot product between s1 and s2.
%---------------------------------------------------------------------------

my_dot=s1(2:Lp,2:Mp).*                                                   ...
       s2(2:Lp,2:Mp).*                                                   ...
       K.rmask(2:Lp,2:Mp);

DotProd=sum(sum(my_dot));

return
