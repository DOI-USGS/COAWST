function [Uout,Vout]=rotate_vec(Uinp,Vinp,angle,irotate)

%
% ROTATE_VEC:  Rotate vector field from (lon,lat) to (XI, ETA) or
%                                       (XI, ETA) to (lon,lat)
%
% [Uout,Vout]=rotate_vec(Uinp,Vinp,angle,irotate);
%
% This function rotates provided vector field from geographical
% coordinates (lon,lat) to curvilinear coordinates (XI,ETA) or
% viceversa.  Vector components in geographical coordinates are
% oriented to TRUE EAST (positive) and TRUE NORTH (positive).
%
% The input vector components (Uinp,Vinp) may be located at either
% the center of the cell (RHO-points) or at the staggered Arakawa's
% C-grid locations.
%
% On Input:
%
%    Uinp        U-component data at RHO- or U-points (array)
%
%    Vinp        V-component data at RHO- or V-points (array)
%
%    angle       Curvilinear rotation angle at RHO-points: angle
%                  between XI-axis and TRUE EAST (radians, 2D array)
%
%    irotate     Vector rotation flag:
%                  
%                  irotate = 0,   rotate from (lon,lat) to (XI, ETA)
%
%                    Uout = Uinp * cos(angle) + Vinp * sin(angle)
%                    Vout = Vinp * cos(angle) - Uinp * sin(angle)
%
%                  irotate = 1,   rotate from (XI, ETA) to (lon,lat)
%
%                    Uout = Uinp * cos(angle) - Vinp * sin(angle)
%                    Vout = Vinp * cos(angle) + Uinp * sin(angle)
%
% On Output:
%
%    Uout        U-component data at RHO- or U-points (array)
%
%    Vout        V-component data at RHO- or V-points (array)
%

% svn $Id: rotate_vec.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Determine if 2D or 3D vector field and dimension size.

[Lp,Mp]=size(angle);

L=Lp-1;  Lm=L-1;  Lm2=Lm-1;
M=Mp-1;  Mm=M-1;  Mm2=Mm-1;

if (length(size(Uinp)) == 2),
  is3d=false;
  [Lu,Mu]=size(Uinp);
  [Lv,Mv]=size(Vinp);
elseif (length(size(Uinp)) == 3),
  is3d=true;
  [Lu,Mu,N]=size(Uinp);
  [Lv,Mv,N]=size(Vinp);
end

%  Determine if input vector components are at RHO-points.

rho_points = false;

if ((Lu == Lp) && (Mv == Mp)),
  rho_points = true;
end

%  Replicate rotation angle to a 3D array.

if (is3d),
  angle=repmat(angle,[1,1,N]);
end

%--------------------------------------------------------------------------
%  Rotate from (lon,lat) to (XI,ETA) coordinates.
%--------------------------------------------------------------------------

if (irotate == 0),

  if (rho_points),
    Uout=Uinp.*cos(angle)+Vinp.*sin(angle);
    Vout=Vinp.*cos(angle)-Uinp.*sin(angle);
  else
    Uout=NaN(size(Uinp));
    Vout=NaN(size(Vinp));
    if (is3d),
      VatU(2:Lp,2:Mm,:)=0.25.*(Vinp(1:L ,1:Mm2,:)+                    ...
                               Vinp(2:Lp,1:Mm2,:)+                    ...
                               Vinp(1:L ,2:Mm ,:)+                    ...
                               Vinp(2:Lp,2:Mm ,:));
      Uangle(2:Lp,2:Mm,:)=0.5.*(angle(1:L ,2:Mm,:)+                   ...
                                angle(2:Lp,2:Mm,:));
      Uout(1:L,2:Mm,:)=Uinp(1:L ,2:Mm,:).*cos(Uangle(2:Lp,2:Mm,:))+   ...
                       VatU(2:Lp,2:Mm,:).*sin(Uangle(2:Lp,2:Mm,:));

      for j=[1,M],
        VatU(2:Lp,1,:)=0.5.*(Vinp (1:L,j,:)+Vinp (2:Lp,j,:));
        Uangle(2:Lp,1,:)=0.5.*(angle(1:L,j,:)+angle(2:Lp,j,:));
        Uout(1:L,j,:)=Uinp(1:L ,j,:).*cos(Uangle(2:Lp,1,:))+          ...
                      VatU(2:Lp,1,:).*sin(Uangle(2:Lp,1,:));
      end

      UatV(2:Lm,2:M,:)=0.25.*(Uinp(1:Lm2,1:Mm,:)+                     ...
                              Uinp(2:Lm ,1:Mm,:)+                     ...
                              Uinp(1:Lm2,2:M ,:)+                     ...
                              Uinp(2:Lm ,2:M ,:));
      Vangle(2:Lm,2:M,:)=0.5*(angle(2:Lm,1:Mm,:)+                     ...
                              angle(2:Lm,2:M ,:));

      Vout(2:Lm,1:Mm,:)=Vinp(2:Lm,1:Mm,:).*cos(Vangle(2:Lm,1:Mm,:))-  ...
                        UatV(2:Lm,2:M ,:).*sin(Vangle(2:Lm,1:Mm,:));
                           
      for i = [1,L],
        UatV(1,2:M,:)=0.5.*(Uinp(i,1:Mm,:)+Uinp(i,2:M,:));
        Vangle(1,2:M,:)=0.5*(angle(i,1:Mm,:)+angle(i,2:M,:));
        Vout(i,1:Mm,:)=-UatV(1,2:M,:).*sin(Vangle(1,1:Mm,:))-         ...
                      Vinp(i,1:Mm,:).*cos(Vangle(1,1:Mm,:));
      end
    
    else

      VatU(2:Lp,2:Mm)=0.25.*(Vinp(1:L ,1:Mm2)+                        ...
                             Vinp(2:Lp,1:Mm2)+                        ...
                             Vinp(1:L ,2:Mm )+                        ...
                             Vinp(2:Lp,2:Mm ));
      Uangle(2:Lp,2:Mm)=0.5.*(angle(1:L ,2:Mm)+                       ...
                              angle(2:Lp,2:Mm));
      Uout(1:L,2:Mm)=Uinp(1:L ,2:Mm).*cos(Uangle(2:Lp,2:Mm))+         ...
                     VatU(2:Lp,2:Mm).*sin(Uangle(2:Lp,2:Mm));
      for j = [1,M],
        VatU(2:Lp,1)=0.5.*(Vinp(1:L,j)+Vinp(2:Lp,j));
        Uangle(2:Lp,1)=0.5.*(angle(1:L,j)+angle(2:Lp,j));
        Uout(1:L,j)=Uinp(1:L ,j).*cos(Uangle(2:Lp,1))+                ...
                    VatU(2:Lp,1).*sin(Uangle(2:Lp,1));
      end

      UatV(2:Lm,2:M)=0.25.*(Uinp(1:Lm2,1:Mm)+                         ...
                            Uinp(2:Lm ,1:Mm)+                         ...
                            Uinp(1:Lm2,2:M )+                         ...
                            Uinp(2:Lm ,2:M ));
      Vangle(2:Lm,2:M)=0.5*(angle(2:Lm,1:Mm)+                         ...
                            angle(2:Lm,2:M ));
      Vout(2:Lm,1:Mm)=Vinp(2:Lm,1:Mm).*cos(Vangle(2:Lm,1:Mm))-        ...
                      UatV(2:Lm,2:M ).*sin(Vangle(2:Lm,1:Mm));
                          
      for i=[1,L],
        UatV(1,2:M)=0.5.*(Uinp(i,1:Mm)+Uinp(i,2:M));
        Vangle(1,2:M)=0.5*(angle(i,1:Mm)+angle(i,2:M));
        Vout(i,1:Mm)=Vinp(i,1:Mm).*cos(Vangle(1,1:Mm))-               ...
                     UatV(1,2:M ).*sin(Vangle(1,1:Mm));
      end
	
    end
  end

%--------------------------------------------------------------------------
%  Rotate from (XI,ETA) to (lon,lat) coordinates.
%--------------------------------------------------------------------------

elseif (irotate == 1),

  if (rho_points),
    Uout=Uinp.*cos(angle)-Vinp.*sin(angle);
    Vout=Vinp.*cos(angle)+Uinp.*sin(angle);
  else
    Uout=NaN(size(Uinp));
    Vout=NaN(size(Vinp));
    if (is3d),
      VatU(2:Lp,2:Mm,:)=0.25.*(Vinp(1:L ,1:Mm2,:)+                    ...
                               Vinp(2:Lp,1:Mm2,:)+                    ...
                               Vinp(1:L ,2:Mm ,:)+                    ...
                               Vinp(2:Lp,2:Mm ,:));
      Uangle(2:Lp,2:Mm,:)=0.5.*(angle(1:L ,2:Mm,:)+                   ...
                                angle(2:Lp,2:Mm,:));
      Uout(1:L,2:Mm,:)=Uinp(1:L ,2:Mm,:).*cos(Uangle(2:Lp,2:Mm,:))-   ...
                       VatU(2:Lp,2:Mm,:).*sin(Uangle(2:Lp,2:Mm,:));

      for j=[1,M],
        VatU(2:Lp,1,:)=0.5.*(Vinp (1:L,j,:)+Vinp (2:Lp,j,:));
        Uangle(2:Lp,1,:)=0.5.*(angle(1:L,j,:)+angle(2:Lp,j,:));
        Uout(1:L,j,:)=Uinp(1:L ,j,:).*cos(Uangle(2:Lp,1,:))-          ...
                      VatU(2:Lp,1,:).*sin(Uangle(2:Lp,1,:));
      end

      UatV(2:Lm,2:M,:)=0.25.*(Uinp(1:Lm2,1:Mm,:)+                     ...
                              Uinp(2:Lm ,1:Mm,:)+                     ...
                              Uinp(1:Lm2,2:M ,:)+                     ...
                              Uinp(2:Lm ,2:M ,:));
      Vangle(2:Lm,2:M,:)=0.5*(angle(2:Lm,1:Mm,:)+                     ...
                              angle(2:Lm,2:M ,:));

      Vout(2:Lm,1:Mm,:)=Vinp(2:Lm,1:Mm,:).*cos(Vangle(2:Lm,1:Mm,:))+  ...
                        UatV(2:Lm,2:M ,:).*sin(Vangle(2:Lm,1:Mm,:));
                           
      for i = [1,L],
        UatV(1,2:M,:)=0.5.*(Uinp(i,1:Mm,:)+Uinp(i,2:M,:));
        Vangle(1,2:M,:)=0.5*(angle(i,1:Mm,:)+angle(i,2:M,:));
        Vout(i,1:Mm)=-UatV(1,2:M ).*sin(Vangle(1,1:Mm))+              ...
                      Vinp(i,1:Mm).*cos(Vangle(1,1:Mm));
      end
    
    else

      VatU(2:Lp,2:Mm)=0.25.*(Vinp(1:L ,1:Mm2)+                        ...
                             Vinp(2:Lp,1:Mm2)+                        ...
                             Vinp(1:L ,2:Mm )+                        ...
                             Vinp(2:Lp,2:Mm ));
      Uangle(2:Lp,2:Mm)=0.5.*(angle(1:L ,2:Mm)+                       ...
                              angle(2:Lp,2:Mm));
      Uout(1:L,2:Mm)=Uinp(1:L ,2:Mm).*cos(Uangle(2:Lp,2:Mm))-         ...
                     VatU(2:Lp,2:Mm).*sin(Uangle(2:Lp,2:Mm));
      for j = [1,M],
        VatU(2:Lp,1)=0.5.*(Vinp(1:L,j)+Vinp(2:Lp,j));
        Uangle(2:Lp,1)=0.5.*(angle(1:L,j)+angle(2:Lp,j));
        Uout(1:L,j)=Uinp(1:L ,j).*cos(Uangle(2:Lp,1))-                ...
                    VatU(2:Lp,1).*sin(Uangle(2:Lp,1));
      end

      UatV(2:Lm,2:M)=0.25.*(Uinp(1:Lm2,1:Mm)+                         ...
                            Uinp(2:Lm ,1:Mm)+                         ...
                            Uinp(1:Lm2,2:M )+                         ...
                            Uinp(2:Lm ,2:M ));
      Vangle(2:Lm,2:M)=0.5*(angle(2:Lm,1:Mm)+                         ...
                            angle(2:Lm,2:M ));
      Vout(2:Lm,1:Mm)=Vinp(2:Lm,1:Mm).*cos(Vangle(2:Lm,1:Mm))+        ...
                      UatV(2:Lm,2:M ).*sin(Vangle(2:Lm,1:Mm));
                          
      for i=[1,L],
        UatV(1,2:M)=0.5.*(Uinp(i,1:Mm)+Uinp(i,2:M));
        Vangle(1,2:M)=0.5*(angle(i,1:Mm)+angle(i,2:M));
        Vout(i,1:Mm)=Vinp(i,1:Mm).*cos(Vangle(1,1:Mm))+               ...
                     UatV(1,2:M ).*sin(Vangle(1,1:Mm));
      end
    end
  
  end

end

return
