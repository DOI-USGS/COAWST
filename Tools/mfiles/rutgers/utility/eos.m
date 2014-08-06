function [F]=eos(gname,fname,tindex);

%
% EOS:  Computes ROMS nonlinear equation of state for seawater fields
%
% [F]=eos(gname,fname,tindex);
%
%
% This function computes quantities associated with ROMS nonlinear
% equation of state for seawater.
%
% On Input:
%
%    gname       ROMS grid NetCDF file name (character string)
%    fname       ROMS history NetCDF file name (character string)
%    tindex      Time index (integer)
%
% On Output:
%
%    F           Equation of state fields (structure array)
%                  F.den      "in situ" density
%                  F.bvf      Brunt-Vaisala frequency
%                  F.alpha    Thermal expansion
%                  F.beta     Saline contraction
%                  F.gamma    Adiabatic and isentropic compressibility
%                  F.svel     Sound speed
%                  F.neutral  Neutral surface coefficient
%                  F.lon      Longitude positions
%                  F.lat      Latitude positions
%                  F.Zr       Depth positions at RHO-points
%                  F.Zw       Depth positions at W-points
%                  F.mask     Land/Sea mask
%
% Check Values: (T=3 C, S=35.5, Z=-5000 m)
%
%     alpha = 2.1014611551470D-04 (1/Celsius)
%     beta  = 7.2575037309946D-04 (nondimensional)
%     gamma = 3.9684764511766D-06 (1/Pa)
%     den   = 1050.3639165364     (kg/m3)
%     den1  = 1028.2845117925     (kg/m3)
%     sound = 1548.8815240223     (m/s)
%     bulk  = 23786.056026320     (Pa)
%
% Note: Salinity does not have physical units. Check the following forum
%       post for details:
%
%       https://www.myroms.org/forum/viewtopic.php?f=14&t=294
%
%  Reference:
%
%  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of
%    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.
%    and Oceanic Techn., vol. 12, pp. 381-389.
%

% svn $Id: eos.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%----------------------------------------------------------------------------
%  Set equation of state expansion coefficients.
%----------------------------------------------------------------------------

A00=+19092.56D0;  A01=+209.8925D0;   A02=-3.041638D0;   A03=-1.852732D-3;
A04=-1.361629D-5; B00=+104.4077D0;   B01=-6.500517D0;   B02=+0.1553190D0;
B03=+2.326469D-4; D00=-5.587545D0;   D01=+0.7390729D0;  D02=-1.909078D-2;
E00=+4.721788D-1; E01=+1.028859D-2;  E02=-2.512549D-4;  E03=-5.939910D-7;
F00=-1.571896D-2; F01=-2.598241D-4;  F02=+7.267926D-6;  G00=+2.042967D-3;
G01=+1.045941D-5; G02=-5.782165D-10; G03=+1.296821D-7;  H00=-2.595994D-7;
H01=-1.248266D-9; H02=-3.508914D-9;  Q00=+999.842594D0; Q01=+6.793952D-2;
Q02=-9.095290D-3; Q03=+1.001685D-4;  Q04=-1.120083D-6;  Q05=+6.536332D-9;
U00=+0.824493D0;  U01=-4.08990D-3;   U02=+7.64380D-5;   U03=-8.24670D-7;
U04=+5.38750D-9;  V00=-5.72466D-3;   V01=+1.02270D-4;   V02=-1.65460D-6;
W00=+4.8314D-4; 

grav=9.81;

%----------------------------------------------------------------------------
%  Read in temperature and salinity.
%----------------------------------------------------------------------------

T=nc_read(fname,'temp',tindex);
S=nc_read(fname,'salt',tindex);

[L,M,N]=size(T);

%----------------------------------------------------------------------------
%  Compute depths at RHO- and W- points.
%----------------------------------------------------------------------------

[Zr]=depths(fname,gname,1,0,tindex);
[Zw]=depths(fname,gname,5,0,tindex);

%----------------------------------------------------------------------------
%  Compute density (kg/m3) at standard one atmosphere pressure.
%----------------------------------------------------------------------------

ind=find(T < -2.0);
if (~isempty(ind)),
  T(ind)=-2.0;                   % lower temperature valid minumum value
end,

ind=find(S < 0.0);
if (~isempty(ind)),
  S(ind)=0.0;                    % lower salinity valid minimum value
end,

clear ind

sqrtS=sqrt(S);

den1 = Q00 + Q01.*T + Q02.*T.^2 + Q03.*T.^3 + Q04.*T.^4 + Q05.*T.^5 + ...
       U00.*S + U01.*S.*T + U02.*S.*T.^2 + U03.*S.*T.^3 + U04.*S.*T.^4 + ...
       V00.*S.*sqrtS + V01.*S.*sqrtS.*T + V02.*S.*sqrtS.*T.^2 + ...
       W00.*S.^2;

%----------------------------------------------------------------------------
%  Compute secant bulk modulus (bulk = K0 - K1*z + K2*z*z).
%----------------------------------------------------------------------------

K0 = A00 + A01.*T + A02.*T.^2 + A03.*T.^3 + A04.*T.^4 + ...
     B00.*S + B01.*S.*T + B02.*S.*T.^2 + B03.*S.*T.^3 + ...
     D00.*S.*sqrtS + D01.*S.*sqrtS.*T + D02.*S.*sqrtS.*T.^2;

K1 = E00 + E01.*T + E02.*T.^2 + E03.*T.^3 + ...
     F00.*S + F01.*S.*T + F02.*S.*T.^2 + ...
     G00.*S.*sqrtS;

K2 = G01 + G02.*T + G03.*T.^2 + ...
     H00.*S + H01.*S.*T + H02.*S.*T.^2;

bulk = K0 - K1.*Zr + K2.*Zr.^2;

%----------------------------------------------------------------------------
%  Compute "in situ" density anomaly (kg/m3).
%----------------------------------------------------------------------------

F.den = (den1.*bulk) ./ (bulk + 0.1.*Zr);

%----------------------------------------------------------------------------
%  Compute thermal expansion (1/Celsius), saline contraction (nondimensional),
%  and adiabatic and isentropic compressibility (1/Pa) coefficients.
%----------------------------------------------------------------------------

%  Compute d(den1)/d(S) and d(den1)/d(T) derivatives.

Dden1DS = U00 + U01.*T + U02.*T.^2 + U03.*T.^3 + U04.*T.^4 + ...
          V00.*1.5.*sqrtS + V01.*1.5.*sqrtS.*T + V02.*1.5.*sqrtS.*T.^2 + ...
          W00.*2.*S;

Dden1DT = Q01 + Q02.*2.*T + Q03.*3.*T.^2 + Q04.*4.*T.^3 + Q05.*5.*T.^4 + ...
          U01.*S + U02.*2.*S.*T + U03.*S.*3.*T.^2 + U04.*S.*4.*T.^3 + ...
          V01.*S.*sqrtS + V02.*S.*sqrtS.*2.*T;

%  Compute d(bulk)/d(S), d(bulk)/d(T), and d(bulk)/d(P) derivatives.

DbulkDS = B00 + B01.*T + B02.*T.^2 + B03.*T.^3 + ...
          D00.*1.5.*sqrtS+ D01.*1.5.*sqrtS.*T + D02.*1.5.*sqrtS.*T.^2 - ...
          F00.*Zr - F01.*Zr.*T - F02.*Zr.*T.^2 - ...
          G00.*Zr.*1.5.*sqrtS + ...
          H00.*Zr.^2 + H01.*Zr.^2.*T + H02.*Zr.^2.*T.^2;


DbulkDT = A01 + A02.*2.*T + A03.*3.*T.^2 + A04.*4*T.^3 + ...
          B01.*S + B02.*S.*2.*T + B03.*S.*3.*T.^2 + ...
          D01.*S.*sqrtS + D02.*S.*sqrtS.*2.*T - ...
          E01.*Zr - E02.*Zr.*2.*T - E03.*Zr.*3.*T.^2 - ...
          F01.*Zr.*S - F02.*Zr.*S.*2.*T + ...
          G02.*Zr.^2 + G03.*Zr.^2.*2.*T + ...
          H01.*Zr.^2.*S + H02.*Zr.^2.*S.*2.*T;

DbulkDP = -K1 + K2.*2.*Zr;

wrk = F.den .* (bulk + 0.1.*Zr).^2;

%  Compute thermal expansion (1/Celsius), saline contraction (nondimensional),
%  and adiabatic and isentropic compressibility (1/Pa) coefficients.

F.alpha = -(DbulkDT.*0.1.*Zr.*den1 + ...
            Dden1DT.*bulk.*(bulk+0.1.*Zr)) ./ wrk;

F.beta  =  (DbulkDS.*0.1.*Zr.*den1 + ...
            Dden1DS.*bulk.*(bulk+0.1.*Zr)) ./ wrk;

F.gamma = -0.1.*den1 .* (DbulkDP.*Zr - bulk) ./ wrk;

clear DbulkDP DbulkDS DbulkDT Dden1DP Dden1DS Dden1DT wrk

%----------------------------------------------------------------------------
%  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
%  and vertical W-points:
%----------------------------------------------------------------------------

for k=1:N-1,

  bulk_up = K0(:,:,k+1) - ...
            Zw(:,:,k+2).*(K1(:,:,k+1)-Zw(:,:,k+2).*K2(:,:,k+1));
  bulk_dn = K0(:,:,k  ) - ...
            Zw(:,:,k+1).*(K1(:,:,k  )-Zw(:,:,k+1).*K2(:,:,k  ));

  den_up = (den1(:,:,k+1).*bulk_up) ./ (bulk_up + 0.1.*Zw(:,:,k+2));
  den_dn = (den1(:,:,k  ).*bulk_dn) ./ (bulk_dn + 0.1.*Zw(:,:,k+1));

  F.bvf(:,:,k) = -grav * (den_up-den_dn) ./ ...
                 (0.5.*(den_up + den_dn) .* (Zr(:,:,k+1)-Zr(:,:,k)));

end,

F.bvf(:,:,N)=F.bvf(:,:,N-1);

clear bulk_up bulk_dn den_up den_dn K0 K1 K2

%----------------------------------------------------------------------------
%  Compute sound speed.
%----------------------------------------------------------------------------

svel2=abs(10000.0./(F.den.*F.gamma));
F.svel=sqrt(svel2);

%----------------------------------------------------------------------------
%  Compute factor to convert isopycnal slopes from "in situ" density
%  to neutral surfaces.
%----------------------------------------------------------------------------

F.neutral=1 + grav.*grav ./ (svel2.*abs(F.bvf));

%----------------------------------------------------------------------------
%  Read in and set horizontal positions.
%----------------------------------------------------------------------------

rmask=nc_read(fname,'mask_rho');
rlon=nc_read(fname,'lon_rho');
rlat=nc_read(fname,'lat_rho');

F.mask=repmat(rmask,[1 1 N]);
F.lon=repmat(rlon,[1 1 N]);
F.lat=repmat(rlat,[1 1 N]);

clear rlon rlat

%----------------------------------------------------------------------------
%  Load depths into structure array.
%----------------------------------------------------------------------------

F.Zr=Zr;
F.Zw=Zw;

return

