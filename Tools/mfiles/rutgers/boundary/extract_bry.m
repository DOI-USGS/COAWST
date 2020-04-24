function [B]=extract_bry(ncfile, vname, rec, compact)

%
% EXTRACT_BRY:  Reads requested variable and extracts boundary edges
%
% [B]=extract_bry(ncfile, vname, rec, flag)
%
% This function reads requested variable from a ROMS NetCDF file at the
% specified time record and extracts the lateral boundary edges.  No
% interpolation is carried out.
%
% On Input:
%
%    ncfile        ROMS NetCDF file name (string)
%    vname         NetCDF variable to process (string)
%    rec           Record to process (scalar)
%    compact       Extraction switch:
%
%                    compact = true,   loads extracted boundary edges
%                                      into single array with an extra
%                                      outer dimension of size 4.
%
%                                      B(:,4)      for 2D fields
%                                      B(:,:,4)    for 3D fields
%
%                                      iwest  = 1  western  edge
%                                      isouth = 2  southern edge
%                                      iwest  = 3  eastern  edge
%                                      iwest  = 4  northern edge
%
%                    compact = false,  loads extracted boundary edges
%                                      into an structure with 4 different
%                                      arrays
%
%                                      B.west      western  edge
%                                      B.south     southern edge
%                                      B.east      eastern  edge
%                                      B.north     northern edge
% On Ouput:
%
%    B             Extracted boundary edges (array or structure)
%
  
% svn $Id: extract_bry.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

iwest  = 1;           % western  edge
isouth = 2;           % southern edge
ieast  = 3;           % eastern  edge
inorth = 4;           % northern edge

% If compact extraction, determine the size of the IorJ dimension.

if (compact)
  D = nc_dinfo(ncfile);

  for n=1:length(D)
    name = char(D(n).Name);
    switch name
      case 'xi_rho'
        Lr = D(n).Length;
      case 'xi_u'
        Lu = D(n).Length;
      case 'xi_v'
        Lv = D(n).Length;
      case 'eta_rho'
        Mr = D(n).Length;
      case 'eta_u'
        Mu = D(n).Length;
      case 'eta_v'
        Mv = D(n).Length;
      case 's_rho'
        Nr = D(n).Length;
      case 's_w'
        Nw = D(n).Length;
    end
  end
  IorJ = max(Lr, Mr);      %  maximum RHO-points value of X- or Y-direction
end  

%----------------------------------------------------------------------------
%  Extract boundary edges into a compact array, used in 4D-Var.
%----------------------------------------------------------------------------

if (compact)

  F = nc_read(ncfile, vname, rec);
  
  if (size(F) == 2)                                            % 2D fields
    [Im, Jm] = size(F);
    B = zeros(IorJ, 4);
    if ((Im == Lu) && (Jm == Mu))                              % U-points
      B(1:Jm  ,iwest ) = squeeze(F(1   ,1:Jm));
      B(2:Im+1,isouth) = squeeze(F(1:Im,1   ));
      B(1:Jm  ,ieast ) = squeeze(F(Im  ,1:Jm));
      B(2:Im+1,inorth) = squeeze(F(1:Im,Jm  ));
    elseif ((Im == Lv) && (Jm == Mv))                          % V-points
      B(2:Jm+1,iwest ) = squeeze(F(1   ,1:Jm));
      B(1:Im  ,isouth) = squeeze(F(1:Im,1   ));
      B(2:Jm+1,ieast ) = squeeze(F(Im  ,1:Jm));
      B(1:Im  ,inorth) = squeeze(F(1:Im,Jm  ));
    else                                                       % RHO-points
      B(1:Jm  ,iwest ) = squeeze(F(1   ,1:Jm));
      B(1:Im  ,isouth) = squeeze(F(1:Im,1   ));
      B(1:Jm  ,ieast ) = squeeze(F(Im  ,1:Jm));
      B(1:Im  ,inorth) = squeeze(F(1:Im,Jm  ));
    end
  else                                                         % 3D fields
    [Im, Jm, Km] = size(F);
    B = zeros(IorJ, Km, 4);
    if ((Im == Lu) && (Jm == Mu))                              % U-points
      B(1:Jm  ,1:Km,iwest ) = squeeze(F(1   ,1:Jm,1:Km));
      B(2:Im+1,1:Km,isouth) = squeeze(F(1:Im,1   ,1:Km));
      B(1:Jm  ,1:Km,ieast ) = squeeze(F(Im  ,1:Jm,1:Km));
      B(2:Im+1,1:Km,inorth) = squeeze(F(1:Im,Jm  ,1:Km));
    elseif ((Im == Lv) && (Jm == Mv))                          % V-points
      B(2:Jm+1,1:Km,iwest ) = squeeze(F(1   ,1:Jm,1:Km));
      B(1:Im  ,1:Km,isouth) = squeeze(F(1:Im,1   ,1:Km));
      B(2:Jm+1,1:Km,ieast ) = squeeze(F(Im  ,1:Jm,1:Km));
      B(1:Im  ,1:Km,inorth) = squeeze(F(1:Im,Jm  ,1:Km));
    else                                                       % RHO-points
      B(1:Jm  ,1:Km,iwest ) = squeeze(F(1   ,1:Jm,1:Km));
      B(1:Im  ,1:Km,isouth) = squeeze(F(1:Im,1   ,1:Km));
      B(1:Jm  ,1:Km,ieast ) = squeeze(F(Im  ,1:Jm,1:Km));
      B(1:Im  ,1:Km,inorth) = squeeze(F(1:Im,Jm  ,1:Km));
    end
  end
    
end

%----------------------------------------------------------------------------
%  Extract boundary edges into a compact array, used in open boundary
%  conditions.
%----------------------------------------------------------------------------

if (~compact)

  F = nc_read(ncfile, vname, rec);

  if (size(F) == 2)                                            % 2D fields
    B.west  = squeeze(F(1  ,:  ));
    B.south = squeeze(F(:  ,1  ));
    B.east  = squeeze(F(end,:  ));
    B.north = squeeze(F(:  ,end));
  else                                                         % 3D fields
    B.west  = squeeze(F(1  ,:  ,:));
    B.south = squeeze(F(:  ,1  ,:));
    B.east  = squeeze(F(end,:  ,:));
    B.north = squeeze(F(:  ,end,:));
  end
end

return
