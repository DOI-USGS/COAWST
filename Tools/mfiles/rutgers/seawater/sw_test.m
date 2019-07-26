function [] = sw_test()

% SW_TEST    Test SEAWATER Library Routines
%=========================================================================
% SW_TEST   $Id: sw_test.m 330 2009-03-10 05:57:42Z arango $
%           Copyright (C) CSIRO, Phil Morgan 1994
%
% sw_test
%
% DESCRIPTION:
%    Execute test routines to test and verify SEAWATER Library routines
%    for your platform.  Prints output to screen and to file sw_test.txt
%
%    Use the "more" command to scroll results to screen
%
% OUTPUT:
%   file sw_test.txt
%
% AUTHOR:  Phil Morgan, Lindsay Pender (Lindsay.Pender@csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%

% Modifications
% 03-12-12. Lindsay Pender, Converted to ITS-90.

delete sw_test.txt
disp('OUTPUT FROM THIS TEST WILL ALSO BE SAVED IN FILE sw_test.txt')
disp(' <enter> to continue...')
pause
reply        = input('Full listing of help for each routine (y/n) ? ','s');
display_help = strcmp(reply,'y') | strcmp(reply,'Y');

format compact
echo off
diary sw_test.txt

disp( '***********************')
disp( '    TEST REPORT    ')
disp( ' ')
disp( ' SEA WATER LIBRARY ')
disp( ' ')
sw_ver
disp( ' ')
disp(['Matlab Version ' version ])
disp( ' ')
disp(['   ' date ''])
disp( '***********************')

disp(' ')

%--------------------------------
% TEST MAIN MODULE  sw_ptmp.m
%      SUB-MODULES  sw_atg.m
%--------------------------------
module     = 'sw_ptmp';
submodules = 'sw_adtg.m';
disp('*************************************')
disp(['**  TESTING MODULE: ' module])
disp(['**  and SUB-MODULE: ' submodules])
disp('*************************************')
if display_help
   eval(['help ' module])
   eval(['help ' submodules])
end %if

% TEST 1 - data from Unesco 1983 p45

T    = [ 0  0  0  0  0  0;
        10 10 10 10 10 10;
    20 20 20 20 20 20;
    30 30 30 30 30 30;
    40 40 40 40 40 40 ]/1.00024;

S    = [25 25 25 35 35 35;
        25 25 25 35 35 35;
    25 25 25 35 35 35 ;
    25 25 25 35 35 35;
    25 25 25 35 35 35 ];

P    = [0 5000 10000 0 5000 10000;
        0 5000 10000 0 5000 10000;
    0 5000 10000 0 5000 10000 ;
    0 5000 10000 0 5000 10000;
    0 5000 10000 0 5000 10000 ];

Pr = [0 0 0 0 0 0];

UN_ptmp =      [ 0  -0.3061  -0.9667   0  -0.3856 -1.0974;
                10   9.3531   8.4684  10   9.2906  8.3643;
            20  19.0438  17.9426  20  18.9985 17.8654;
            30  28.7512  27.4353  30  28.7231 27.3851;
            40  38.4607  36.9254  40  38.4498 36.9023];

ptmp    = sw_ptmp(S,T,P,Pr)*1.00024;

%----------------
% DISPLAY RESULTS
%----------------
disp(' ')
disp   ('********************************************************')
disp   ('Comparison of accepted values from UNESCO 1983 ')
disp   (' (Unesco Tech. Paper in Marine Sci. No. 44, p45)')
disp   (['with computed results from ' module ' on ' computer ' computer'])
disp   ('********************************************************')

for icol = 1:length(S(1,:))
disp(' ')
disp   ('   Sal  Temp  Press     PTMP       sw_ptmp')
disp   ('  (psu)  (C)   (db)     (C)          (C)')
  fprintf(1,' %4.0f  %4.0f   %5.0f   %8.4f  %11.5f\n', ...
  [S(:,icol) T(:,icol) P(:,icol) UN_ptmp(:,icol) ptmp(:,icol)]');
end %for

%-------------------------------------------------------------------------------
% TEST MAIN MODULE  sw_svan.m
%      SUB-MODULES  sw_dens.m sw_dens0.m sw_smow.m sw_seck.m sw_pden.m sw_ptmp.m
%------------------------------------------------------------------------------
module     = 'sw_svan.m';
submodules = 'sw_dens.m sw_dens0.m sw_smow.m sw_seck.m sw_pden.m sw_ptmp.m';
disp(' ')
disp('************************************************************************')
disp(['**  TESTING MODULE: ' module])
disp(['**  and SUB-MODULE: ' submodules])
disp('************************************************************************')
if display_help
   eval(['help ' module])
   eval(['help ' submodules])
end %if

% TEST DATA FROM
% Unesco Tech. Paper in Marine Sci. No. 44, p22

s = [0     0   0      0  35    35  35    35]';
p = [0 10000   0  10000   0 10000   0 10000]';
t = ([0     0  30     30   0     0  30    30] / 1.00024)';

UN_svan = [2749.54 2288.61 3170.58 3147.85 ...
              0.0     0.00  607.14  916.34]';

svan    = sw_svan(s,t,p);

%----------------
% DISPLAY RESULTS
%----------------
disp(' ')
disp   ('********************************************************')
disp   ('Comparison of accepted values from UNESCO 1983')
disp   (' (Unesco Tech. Paper in Marine Sci. No. 44, p22)')
disp   (['with computed results from ' module ' on ' computer ' computer'])
disp   ('********************************************************')
disp(' ')
disp   ('   Sal  Temp  Press        SVAN        sw_svan')
disp   ('  (psu)  (C)   (db)    (1e-8*m3/kg)  (1e-8*m3/kg)')
fprintf(1,' %4.0f  %4.0f   %5.0f   %11.2f    %11.3f\n',[s t p UN_svan 1e+8*svan]');

%-------------------------------------------------------------------------------
% TEST MAIN MODULE
%      SUB-MODULES
%------------------------------------------------------------------------------
module     = 'sw_salt.m';
submodules = 'sw_salrt.m sw_salrp.m sw_sals.m';
disp(' ')
disp('************************************************************************')
disp(['**  TESTING MODULE: ' module])
disp(['**  and SUB-MODULE: ' submodules])
disp('************************************************************************')
if display_help
   eval(['help ' module])
   eval(['help ' submodules])
end %if

% TEST 1 - data from Unesco 1983 p9
%***************************************************************************

R     = [1 1.2 0.65]';   % cndr = R
T     = ([15 20 5]/1.00024)';
P     = [0 2000 1500]';

Rt    = [1 1.0568875 0.81705885]';
UN_S  = [35 37.245628 27.995347]';
S     = sw_salt(R,T,P);

%----------------
% DISPLAY RESULTS
%----------------
disp(' ')
disp   ('********************************************************')
disp   ('Comparison of accepted values from UNESCO 1983 ')
disp   (' (Unesco Tech. Paper in Marine Sci. No. 44, p9)')
disp   (['with computed results from ' module ' on ' computer ' computer'])
disp   ('********************************************************')
disp(' ')
disp   ('   Temp    Press       R              S           sw_salt')
disp   ('   (C)     (db)    (no units)       (psu)          (psu) ')
table = [T P R UN_S S]';
fprintf(1,' %4.0f       %4.0f  %8.2f      %11.6f  %14.7f\n', table);

%-------------------------------------------------------------------------------
% TEST MAIN MODULE
%      SUB-MODULES
%------------------------------------------------------------------------------
module     = 'sw_cndr.m';
submodules = 'sw_salds.m';
disp(' ')
disp('************************************************************************')
disp(['**  TESTING MODULE: ' module])
disp(['**  and SUB-MODULE: ' submodules])
disp('************************************************************************')
if display_help
   eval(['help ' module])
   eval(['help ' submodules])
end %if

% TEST 1 - data from Unesco 1983 p9

T    = ([0   10     0   10  10  30]/1.00024)';
P    = [0    0  1000 1000   0   0]';
S    = [25  25    25   25  40  40]';
UN_R = [ 0.498088 0.654990 0.506244 0.662975 1.000073 1.529967]';
R    = sw_cndr(S,T,P);

%----------------
% DISPLAY RESULTS
%----------------
disp(' ')
disp   ('********************************************************')
disp   ('Comparison of accepted values from UNESCO 1983 ')
disp   (' (Unesco Tech. Paper in Marine Sci. No. 44, p14)')
disp   (['with computed results from ' module ' on ' computer ' computer'])
disp   ('********************************************************')
disp(' ')
disp   ('   Temp    Press       S            cndr         sw_cndr')
disp   ('   (C)     (db)      (psu)        (no units)    (no units) ')
table = [T P S UN_R R]';
fprintf(1,' %4.0f       %4.0f   %8.6f   %11.6f  %14.8f\n', table);

%-------------------------------------------------------------------------------
% TEST MAIN MODULE
%      SUB-MODULES
%------------------------------------------------------------------------------
module     = 'sw_dpth.m';
disp(' ')
disp('************************************************************************')
disp(['**  TESTING MODULE: ' module])
disp('************************************************************************')
if display_help
   eval(['help ' module])
end %if

% TEST DATA - matrix "pressure", vector "lat"  Unesco 1983 data p30.

lat = [0 30 45 90];
P   = [  500   500   500   500;
        5000  5000  5000  5000;
       10000 10000 10000 10000];

UN_dpth = [   496.65   496.00   495.34   494.03;
             4915.04  4908.56  4902.08  4889.13;
         9725.47  9712.65  9699.84  9674.23];

dpth = sw_dpth(P,lat);

%----------------
% DISPLAY RESULTS
%----------------
disp(' ')
disp   ('********************************************************')
disp   ('Comparison of accepted values from Unesco 1983 ')
disp   (' (Unesco Tech. Paper in Marine Sci. No. 44, p28)')
disp   (['with computed results from ' module ' on ' computer ' computer'])
disp   ('********************************************************')

for irow = 1:3
   disp(' ')
   disp   ('    Lat       Press     DPTH      sw_dpth')
   disp   ('  (degree)    (db)     (meter)    (meter)')
   table = [lat' P(irow,:)' UN_dpth(irow,:)' dpth(irow,:)'];
   fprintf(1,'  %6.3f     %6.0f   %8.2f   %8.3f\n', table')
end %for

%-------------------------------------------------------------------------------
% TEST MAIN MODULE
%      SUB-MODULES
%------------------------------------------------------------------------------
module     = 'sw_fp.m';
disp(' ')
disp('************************************************************************')
disp(['**  TESTING MODULE: ' module])
disp('************************************************************************')
if display_help
   eval(['help ' module])
end %if

% TEST 1 -
% UNESCO DATA p.30
%***************************************************************************
S    = [ 5   10  15  20  25  30  35  40;
         5   10  15  20  25  30  35  40];

P    = [  0   0   0   0   0  0    0   0;
        500 500 500 500 500 500 500 500];


UN_fp = [-0.274 -0.542 -0.812 -1.083 -1.358 -1.638 -1.922 -2.212;
         -0.650 -0.919 -1.188 -1.460 -1.735 -2.014 -2.299 -2.589];

fp    = sw_fp(S,P);

%----------------
% DISPLAY RESULTS
%----------------
disp(' ')
disp   ('********************************************************')
disp   ('Comparison of accepted values from UNESCO 1983 ')
disp   (' (Unesco Tech. Paper in Marine Sci. No. 44, p30)')
disp   (['with computed results from ' module ' on ' computer ' computer'])
disp   ('********************************************************')

for irow = 1:2
  disp(' ')
  disp   ('   Sal   Press      fp        sw_fp')
  disp   ('  (psu)   (db)      (C)        (C)')
  table = [S(irow,:); P(irow,:); UN_fp(irow,:); fp(irow,:)];
  fprintf(1,' %4.0f   %5.0f   %8.3f  %11.4f\n', table)
end %for

%-------------------------------------------------------------------------------
% TEST MAIN MODULE
%      SUB-MODULES
%------------------------------------------------------------------------------
module     = 'sw_cp.m';
disp(' ')
disp('************************************************************************')
disp(['**  TESTING MODULE: ' module])
disp('************************************************************************')
if display_help
   eval(['help ' module])
end %if

% TEST 1 -
% DATA FROM POND AND PICKARD INTRO. DYNAMICAL OCEANOGRAPHY 2ND ED. 1986
%***************************************************************************

T    = [ 0  0  0  0  0  0;
        10 10 10 10 10 10;
    20 20 20 20 20 20;
    30 30 30 30 30 30;
    40 40 40 40 40 40 ] / 1.00024;

S    = [25 25 25 35 35 35;
        25 25 25 35 35 35;
    25 25 25 35 35 35 ;
    25 25 25 35 35 35;
    25 25 25 35 35 35 ];

P    = [0 5000 10000 0 5000 10000;
        0 5000 10000 0 5000 10000;
    0 5000 10000 0 5000 10000 ;
    0 5000 10000 0 5000 10000;
    0 5000 10000 0 5000 10000 ];

UN_cp =      [  4048.4  3896.3  3807.7  3986.5  3849.3  3769.1;
                4041.8  3919.6  3842.3  3986.3  3874.7  3804.4;
            4044.8  3938.6  3866.7  3993.9  3895.0  3828.3;
            4049.1  3952.0  3883.0  4000.7  3909.2  3844.3;
            4051.2  3966.1  3905.9  4003.5  3923.9  3868.3 ];

cp    = sw_cp(S,T,P);

%----------------
% DISPLAY RESULTS
%----------------
disp(' ')
disp   ('********************************************************')
disp   ('Comparison of accepted values from UNESCO 1983 ')
disp   (' (Unesco Tech. Paper in Marine Sci. No. 44, p37)')
disp   (['with computed results from ' module ' on ' computer ' computer'])
disp   ('********************************************************')

for icol = 1:length(S(1,:))
  disp(' ')
  disp   ('   Sal  Temp  Press      Cp        sw_cp')
  disp   ('  (psu)  (C)   (db)    (J/kg.C)   (J/kg.C)')
  fprintf(1,' %4.0f  %4.0f   %5.0f   %8.1f  %11.2f\n', ...
  [S(:,icol) T(:,icol) P(:,icol) UN_cp(:,icol) cp(:,icol)]');
end %for

%-------------------------------------------------------------------------------
% TEST MAIN MODULE
%      SUB-MODULES
%------------------------------------------------------------------------------
module     = 'sw_svel.m';
disp(' ')
disp('************************************************************************')
disp(['**  TESTING MODULE: ' module])
disp('************************************************************************')
if display_help
   eval(['help ' module])
end %if

% TEST 1 -
% DATA FROM POND AND PICKARD INTRO. DYNAMICAL OCEANOGRAPHY 2ND ED. 1986
%***************************************************************************

T    = [ 0  0  0  0  0  0;
        10 10 10 10 10 10;
    20 20 20 20 20 20;
    30 30 30 30 30 30;
    40 40 40 40 40 40 ] / 1.00024;

S    = [25 25 25 35 35 35;
        25 25 25 35 35 35;
    25 25 25 35 35 35 ;
    25 25 25 35 35 35;
    25 25 25 35 35 35 ];

P    = [0 5000 10000 0 5000 10000;
        0 5000 10000 0 5000 10000;
    0 5000 10000 0 5000 10000 ;
    0 5000 10000 0 5000 10000;
    0 5000 10000 0 5000 10000 ];

UN_svel =      [1435.8  1520.4  1610.4  1449.1  1534.0  1623.2;
                1477.7  1561.3  1647.4  1489.8  1573.4  1659.0;
            1510.3  1593.6  1676.8  1521.5  1604.5  1687.2;
            1535.2  1619.0  1700.6  1545.6  1629.0  1710.1;
            1553.4  1638.0  1719.2  1563.2  1647.3  1727.8 ];

svel    = sw_svel(S,T,P);

%----------------
% DISPLAY RESULTS
%----------------
disp(' ')
disp   ('********************************************************')
disp   ('Comparison of accepted values from UNESCO 1983 ')
disp   (' (Unesco Tech. Paper in Marine Sci. No. 44, p50)')
disp   (['with computed results from ' module ' on ' computer ' computer'])
disp   ('********************************************************')

for icol = 1:length(S(1,:))
  disp(' ')
  disp   ('   Sal  Temp  Press     SVEL       sw_svel')
  disp   ('  (psu)  (C)   (db)     (m/s)       (m/s)')
  fprintf(1,' %4.0f  %4.0f   %5.0f   %8.1f  %11.3f\n', ...
  [S(:,icol) T(:,icol) P(:,icol) UN_svel(:,icol) svel(:,icol)]');
end %for

%----------------------------------------------------------------------------
% TEST MAIN MODULE
%      SUB-MODULES
%---------------------------------------------------------------------------
submodules     = 'sw_alpha.m sw_beta.m sw_aonb.m';
disp(' ')
disp('**********************************************************************')
disp(['**  TESTING MODULE: ' submodules])
disp('**********************************************************************')
if display_help
   eval(['help ' submodules])
end %if

% DATA FROM MCDOUOGALL 1987
s    = 40;
ptmp = 10;
p    = 4000;
beta_lit  = 0.72088e-03;
aonb_lit  = 0.34763;
alpha_lit = aonb_lit*beta_lit;

%$$$ % TEST ARGUMENT PASSING
%$$$ beta = sw_beta(s,ptmp,p,'ptmp')
%$$$ beta = sw_beta(s,ptmp,p,'temp')
%$$$ beta = sw_beta(s,ptmp,p)
%$$$
%$$$ alpha = sw_alpha(s,ptmp,p,'ptmp')
%$$$ alpha = sw_alpha(s,ptmp,p,'temp')
%$$$ alpha = sw_alpha(s,ptmp,p)
%$$$
%$$$ aonb  = sw_aonb( s,ptmp,p,'ptmp')
%$$$ aonb  = sw_aonb( s,ptmp,p,'temp')
%$$$ aonb  = sw_aonb( s,ptmp,p)
%$$$
beta  = sw_beta( s,ptmp,p,'ptmp');
alpha = sw_alpha(s,ptmp,p,'ptmp');
aonb  = sw_aonb( s,ptmp,p,'ptmp');

%----------------
% DISPLAY RESULTS
%----------------
disp(' ')
disp   ('********************************************************')
disp   ('Comparison of accepted values from MCDOUGALL 1987 ')
disp   (['with computed results on ' computer ' computer'])
disp   ('********************************************************')

  disp(' ')
  disp   ('   Sal  Temp  Press     BETA       sw_beta')
  disp   ('  (psu)  (C)   (db)   (psu^-1)     (psu^-1)')
  fprintf(1,' %4.0f  %4.0f   %5.0f   %11.4e  %11.5e\n', ...
  [s ptmp p beta_lit beta]');

  disp(' ')
  disp   ('   Sal  Temp  Press     AONB       sw_aonb')
  disp   ('  (psu)  (C)   (db)   (psu C^-1)   (psu C^-1)')
  fprintf(1,' %4.0f  %4.0f   %5.0f   %8.5f  %11.6f\n', ...
  [s ptmp p aonb_lit aonb]');

  disp(' ')
  disp   ('   Sal  Temp  Press     ALPHA       sw_alpha')
  disp   ('  (psu)  (C)   (db)    (psu^-1)     (psu^-1)')
  fprintf(1,' %4.0f  %4.0f   %5.0f   %11.4e  %11.4e\n', ...
  [s ptmp p alpha_lit alpha]');

%--------------------------------
% TEST MAIN MODULE  sw_satO2.m
%      SUB-MODULES
%--------------------------------
module     = 'sw_satO2 sw_satN2 sw_satAr';
disp(' ')
disp('*************************************')
disp(['**  TESTING MODULE: ' module])
disp(['**  and SUB-MODULE: ' submodules])
disp('*************************************')
if display_help
   eval(['help ' module])
end %if

% Data from Weiss 1970

T    = [-1 -1;
        10 10;
    20 20 ;
    40 40 ] / 1.00024;

S    = [20 40;
        20 40;
    20 40 ;
    20 40];

lit_O2=  [ 9.162   7.984;
           6.950   6.121;
       5.644   5.015;
       4.050   3.656];

lit_N2=  [16.28   14.01;
          12.64   11.01;
      10.47    9.21;
       7.78    6.95];

lit_Ar=  [ 0.4456 0.3877;
           0.3397 0.2989;
       0.2766 0.2457;
       0.1986 0.1794];


satO2    = sw_satO2(S,T);
satN2    = sw_satN2(S,T);
satAr    = sw_satAr(S,T);

%----------------
% DISPLAY RESULTS
%----------------
disp(' ')
disp   ('********************************************************')
disp   ('Comparison of accepted values from Weiss, R.F. 1979 ')
disp   ('"The solubility of nitrogen, oxygen and argon in water and seawater."')
disp   (' Deap-Sea Research., 1970, Vol 17, pp721-735.')
disp   (['with computed results from ' module ' on ' computer ' computer'])
disp   ('********************************************************')

for icol = 1:length(S(1,:))
disp(' ')
disp   ('   Sal  Temp      O2         sw_satO2')
disp   ('  (psu)  (C)      (ml/l)     (ml/l)')
  fprintf(1,' %4.0f  %4.0f    %8.2f   %9.3f\n', ...
  [S(:,icol) T(:,icol)  lit_O2(:,icol) satO2(:,icol)]');
end %for

for icol = 1:length(S(1,:))
disp(' ')
disp   ('   Sal  Temp      N2         sw_satN2')
disp   ('  (psu)  (C)      (ml/l)     (ml/l)')
  fprintf(1,' %4.0f  %4.0f    %8.2f  %9.3f\n', ...
  [S(:,icol) T(:,icol)  lit_N2(:,icol) satN2(:,icol)]');
end %for

for icol = 1:length(S(1,:))
disp(' ')
disp   ('   Sal  Temp      Ar         sw_satAr')
disp   ('  (psu)  (C)      (ml/l)     (ml/l)')
  fprintf(1,' %4.0f  %4.0f     %8.4f  %9.4f\n', ...
  [S(:,icol) T(:,icol)  lit_Ar(:,icol) satAr(:,icol)]');
end %for

diary off

