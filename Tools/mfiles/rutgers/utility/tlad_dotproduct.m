function [dot0, dot1] = tlad_dotproduct(TLname, ADname)

%
% TLAD_DOTPRODUCT: Computes initial and final times TL/AD dot products.
%
% [dot0, dot1] = tlad_dotproduct(TLname, ADname)
%
% This function is used to check the symmetry of ROMS tangent linear (TL)
% and adjoint (AD) operators by computing the initial and final dot product
% of the TL and AD state vectors.
%
% Regardless of the same initial state vectors for the TL ROMS at t=Ti and
% the AD ROMS at t=Tf, the initial and final TL/AD dot products ate the
% same within roundoff.
%
% Let x(t) denote the TL state vector and z(t) denote the adjoint state
% vector solutions. Therefore,
%
%    x(t) = M x(0)	        (1)    TL advanced forward   from t=Ti:Tf
%
%    z(0) = transpose(M) z(t)	(2)    AD advanced backwards from t=Tf:Ti
%
% It is straight forward to prove that:
%
%    traspose[z(t)] * x(t) = tranpose[x(0)] * z(0)
%                     dot1 = dot0
%
% On Input:
%
%    TLname    TL ROMS history NetCDF filename (string)
%    ADname    AD ROMS history NetCDF filename (string)
%
% On Output:
%
%    dot0      Initial TL/AD dot product at time t=Ti
%    dot1      Final   TL/AD dot product at time t=Tf
%
% Here, both the TL ROMS and AD ROMS are initialized with the same state
% vector (usually, random values):
%
%    x(t=Ti) = z(t=Tf)
%
% ROMS needs to be configure that write the history solutions at every
% timestep or at less at the initial and final times.

% svn $Id: tlad_dotproduct.m 1156 2023-02-18 01:44:37Z arango $
%===========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license           Hernan G. Arango        %
%    See License_ROMS.txt                           Andrew M. Moore         %
%===========================================================================%  
  
% Read solution time and set number of TL/AD records.

tl_time = nc_read(TLname, 'ocean_time');
NrecTL = length(tl_time);

ad_time = nc_read(ADname, 'ocean_time');
NrecAD = length(ad_time);

TimeAtt = nc_getatt(TLname, 'units', 'ocean_time');
ind = findstr(TimeAtt,'since');
epoch = datenum(TimeAtt(ind+5:end));

% Set state vector variables to process.

S = nc_vnames(TLname);
vars = {'zeta', 'u', 'v', 'temp'};
if (any(strcmp({S.Variables.Name}, 'salt')))
  vars = [vars, 'salt'];
end

% Use NaN to remove _FillValue due to the land/mask.

disp(blanks(2));

for tlrec = 1:NrecTL
  adrec = NrecTL - tlrec + 1;
  
  tldate = datestr(epoch+tl_time(tlrec)/86400, 'mm-dd HH:MM:SS');
  addate = datestr(epoch+tl_time(adrec)/86400, 'mm-dd HH:MM:SS');
  
  dot = 0;

  for var = vars
    field = char(var);

    tl_f = nc_read(TLname, field, tlrec, NaN);
    ad_f = nc_read(ADname, field, adrec, NaN);
  
    ind = ~isnan(tl_f);

    dot = dot + sum(transpose(tl_f(ind)) * ad_f(ind));
  end

  disp(['TL/AD date: ', tldate, ' / ', addate,                          ...
        ', TLrec = ', num2str(tlrec, '%2.2i'),                          ...
        ', ADrec = ', num2str(adrec, '%2.2i'),                          ...
        ', Dot Product = ', num2str(dot, '%15.8f')]);

  if (tlrec == 1)
    dot0 = dot;
  elseif (tlrec == NrecTL)
    dot1 = dot;
  end

end

disp(blanks(2));
disp(['Initial, dot0 = ', num2str(dot0, '%15.8f')]);
disp(['Final,   dot1 = ', num2str(dot1, '%15.8f')]);
disp(['    dot1-dot0 = ', num2str(dot1-dot0, '%15.8e')]);
