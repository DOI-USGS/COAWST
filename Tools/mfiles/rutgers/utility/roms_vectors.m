function [Uout,Vout]=roms_vectors(Uinp,Vinp,angle,umask,vmask,boundary);

%
% ROMS_VECTORS:  Proccess ROMS vectors from input RHO-point data
%
% [Uout,Vout]=roms_vectors(Uinp,Vinp,angle,umask,vmask,boundary)
%
% This function processes ROMS vector data for either the full
% grid or boundary edges. This function must be used when the
% input vector data (Uinp,Vinp) is at ROMS points.
%
% The strategy is to get any horizontal vector field (Uinp,Vinp)
% at RHO-points for the event that a rotation to ROMS curvilinear
% grid is needed.  If required, the output vectors are computed
% at the appropriate Arakawa C-grid location.
%
% On Input:
%
%    Uinp        U-component data at RHO-points (array or structure)
%    Vinp        V-component data at RHO-points (array or structure)
%    angle       ROMS curvilinear rotation angle at RHO-points
%                  (radians, 2D array)
%    umask       Land/Sea mask at U-points (2D array)
%    vmask       Land/Sea mask at V-points (2D array)
%
%    boundary    ROMS open boundary switch (1D array, OPTIONAL)
%
%                  boundary(1)        Process western  boundary
%                  boundary(2)        Process eastern  boundary
%                  boundary(3)        Process southern boundary
%                  boundary(4)        Process northern boundary
%
% On Output:
%
%    Uout        U-component data at U-points (2D or 3D array)
%    Vout        V-component data at V-points (2D or 3D array)
%
% If processing boundary data, the velocity data are structures
% with the following fields:
%
%    Uinp.west   Vinp.west   Uout.west   Vout.west
%    Uinp.east   Vinp.east   Uout.east   Vout.east
%    Uinp.south  Vinp.south  Uout.south  Vout.south
%    Uinp.north  Vinp.norht  Uout.north  Vout.north
%
% Curvilinear grid rotation to (XI,ETA) coordinates: (RHO-points)
%
%    Urot(i,j,:) =   Uinp(i,j,:) * cos(angle) + Vinp(i,j,:) * sin(angle)
%    Vrot(i,j,:) = - Uinp(i,j,:) * sin(angle) + Vinp(i,j,:) * cos(angle)
%
%
% ROMS staggered, horizontal grid stencil:
%
%
%          --------v(i,j+1,k)---------
%          |                         |
%       u(i,j,k)    r(i,j,k)     u(i+1,j,k)
%          |                         |
%          ---------v(i,j,k)----------
%
%    Uout(i,j,k) = 0.5 * (Uinp(i-1,j,k) + Uinp(i,j,k))
%    Vout(i,j,k) = 0.5 * (Vinp(i,j-1,k) + Vinp(i,j,k))
%

% svn $Id: roms_vectors.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set internal parameters.

if (nargin < 6),
  DoBoundary=0;
else
  DoBoundary=1;
end

if (min(min(angle)) == 0 & max(max(angle)) == 0),
  rotate=0;
else
  rotate=1;
end

if (DoBoundary),
  names=fieldnames(Uinp);
  is2d=0;
  if (length(size(getfield(Uinp,char(names(1))))) == 2),
    is2d=1;
  end

  Bnames={'west','east','south','north'};
  ufields=isfield(Uinp,Bnames);
  vfields=isfield(Vinp,Bnames);
  for ib=1:4,
    if (boundary(ib)),
      if (~ufields(ib)),
        error(['ROMS_VECTORS - cannot find field: Uinp.',             ...
               char(Bnames(ib))]);
      end
      if (~vfields(ib)),
        error(['ROMS_VECTORS - cannot find field: Vinp.',             ...
               char(Bnames(ib))]);
      end
    end
  end
else
  if (length(size(Uinp)) == 2),
    is2d=1;
    [Lr,Mr]=size(Uinp);
  else
    is2d=0;
    [Lr,Mr,N]=size(Uinp);
  end
  Lu=Lr-1;  Mu=Mr;
  Lv=Lr;    Mv=Mr-1;
end

%--------------------------------------------------------------------------
%  Rotate and average to staggered C-grid locations.
%--------------------------------------------------------------------------

if (rotate),
  if (DoBoundary),
    [Lr,Mr]=size(angle);
    for ib=1:4,
      if (boundary(ib)),
        switch ib
          case {1}
            bangle=angle(1:2,:);
            if (is2d),
              Urot= Uinp.west.*cos(bangle)+Vinp.west.*sin(bangle);
              Vrot=-Uinp.west.*sin(bangle)+Vinp.west.*cos(bangle);

              Uout.west=squeeze(0.5.*(Urot(1,:)+Urot(2,:)));
              Vout.west=squeeze(0.5.*(Vrot(1,1:Mr-1)+Vrot(1,2:Mr)));

              Uout.west=squeeze(umask(1,:)').*Uout.west;
              Vout.west=squeeze(vmask(1,:)').*Vout.west;
            else
              N=size(Uinp.west,3);
              bangle=repmat(bangle,[1,1,N]);

              Urot= Uinp.west.*cos(bangle)+Vinp.west.*sin(bangle);
              Vrot=-Uinp.west.*sin(bangle)+Vinp.west.*cos(bangle);

              Uout.west=squeeze(0.5.*(Urot(1,:,:)+Urot(2,:,:)));
              Vout.west=squeeze(0.5.*(Vrot(1,1:Mr-1,:)+Vrot(1,2:Mr,:)));
            
              Uout.west=repmat(squeeze(umask(1,:)'),[1,N]).*Uout.west;
              Vout.west=repmat(squeeze(vmask(1,:)'),[1,N]).*Vout.west;
            end,
          case {2}
            bangle=angle(end-1:end,:);
            if (is2d),
              Urot= Uinp.east.*cos(bangle)+Vinp.east.*sin(bangle);
              Vrot=-Uinp.east.*sin(bangle)+Vinp.east.*cos(bangle);

              Uout.east=squeeze(0.5.*(Urot(1,:)+Urot(2,:)));
              Vout.east=squeeze(0.5.*(Vrot(1,1:Mr-1)+Vrot(1,2:Mr)));

              Uout.east=squeeze(umask(end,:)').*Uout.east;
              Vout.east=squeeze(vmask(end,:)').*Vout.east;
            else
              N=size(Uinp.east,3);
              bangle=repmat(bangle,[1,1,N]);

              Urot= Uinp.east.*cos(bangle)+Vinp.east.*sin(bangle);
              Vrot=-Uinp.east.*sin(bangle)+Vinp.east.*cos(bangle);

              Uout.east=squeeze(0.5.*(Urot(1,:,:)+Urot(2,:,:)));
              Vout.east=squeeze(0.5.*(Vrot(2,1:Mr-1,:)+Vrot(2,2:Mr,:)));

              Uout.east=repmat(squeeze(umask(end,:)'),[1,N]).*Uout.east;
              Vout.east=repmat(squeeze(vmask(end,:)'),[1,N]).*Vout.east;
            end
          case {3}
            bangle=angle(:,1:2);
            if (is2d),
              Urot= Uinp.south.*cos(bangle)+Vinp.south.*sin(bangle);
              Vrot=-Uinp.south.*sin(bangle)+Vinp.south.*cos(bangle);

              Uout.south=squeeze(0.5.*(Urot(1:Lr-1,1)+Urot(2:Lr,1)));
              Vout.south=squeeze(0.5.*(Vrot(:,1)+Vrot(:,2)));

              Uout.south=squeeze(umask(:,1)).*Uout.south;
              Vout.south=squeeze(vmask(:,1)).*Vout.south;
            else
              N=size(Uinp.south,3);
              bangle=repmat(bangle,[1,1,N]);

              Urot= Uinp.south.*cos(bangle)+Vinp.south.*sin(bangle);
              Vrot=-Uinp.south.*sin(bangle)+Vinp.south.*cos(bangle);

              Uout.south=squeeze(0.5.*(Urot(1:Lr-1,1,:)+Urot(2:Lr,1,:)));
              Vout.south=squeeze(0.5.*(Vrot(:,1,:)+Vrot(:,2,:)));

              Uout.south=repmat(squeeze(umask(:,1)),[1,N]).*Uout.south;
              Vout.south=repmat(squeeze(vmask(:,1)),[1,N]).*Vout.south;
            end
          case {4}
            bangle=angle(:,end-1:end);
            if (is2d),
              Urot= Uinp.north.*cos(bangle)+Vinp.north.*sin(bangle);
              Vrot=-Uinp.north.*sin(bangle)+Vinp.north.*cos(bangle);

              Uout.north=squeeze(0.5.*(Urot(1:Lr-1,2)+Urot(2:Lr,2)));
              Vout.north=squeeze(0.5.*(Vrot(:,1)+Vrot(:,2)));

              Uout.north=squeeze(umask(:,end)).*Uout.north;
              Vout.north=squeeze(vmask(:,end)).*Vout.north;
            else
              N=size(Uinp.north,3);
              bangle=repmat(bangle,[1,1,N]);

              Urot= Uinp.north.*cos(bangle)+Vinp.north.*sin(bangle);
              Vrot=-Uinp.north.*sin(bangle)+Vinp.north.*cos(bangle);

              Uout.north=squeeze(0.5.*(Urot(1:Lr-1,2,:)+Urot(2:Lr,2,:)));
              Vout.north=squeeze(0.5.*(Vrot(:,1,:)+Vrot(:,2,:)));
            
              Uout.north=repmat(squeeze(umask(:,end)),[1,N]).*Uout.north;
              Vout.north=repmat(squeeze(vmask(:,end)),[1,N]).*Vout.north;
            end
        end
      end
    end
  else
    if (is2d),
      Urot= Uinp.*cos(angle)+Vinp.*sin(angle);
      Vrot=-Uinp.*sin(angle)+Vinp.*cos(angle);

      Uout=0.5.*(Urot(1:Lr-1,1:Mr)+Urot(2:Lr,1:Mr));
      Vout=0.5.*(Vrot(1:Lr,1:Mr-1)+Vrot(1:Lr,2:Mr));

      Uout=umask.*Uout;
      Vout=vmask.*Vout;
    else
      angle=repmat(angle,[1,1,N]);

      Urot= Uinp.*cos(angle)+Vinp.*sin(angle);
      Vrot=-Uinp.*sin(angle)+Vinp.*cos(angle);

      Uout=0.5.*(Urot(1:Lr-1,1:Mr,:)+Urot(2:Lr,1:Mr,:));
      Vout=0.5.*(Vrot(1:Lr,1:Mr-1,:)+Vrot(1:Lr,2:Mr,:));

      Uout=repmat(umask,[1,1,N]).*Uout;
      Vout=repmat(vmask,[1,1,N]).*Vout;
    end
  end
end

%--------------------------------------------------------------------------
%  Average to staggered C-grid locations, no rotation.
%--------------------------------------------------------------------------

if (~rotate),
  if (DoBoundary),
    [Lr,Mr]=size(angle);
    for ib=1:4,
      if (boundary(ib)),
        switch ib
          case {1}
            if (is2d),
              Uout.west=squeeze(0.5.*(Uinp.west(1,:)+                 ...
                                      Uinp.west(2,:)));
              Vout.west=squeeze(0.5.*(Vinp.west(1,1:Mr-1)+            ...
                                      Vinp.west(1,2:Mr  )));

              Uout.west=squeeze(umask(1,:)').*Uout.west;
              Vout.west=squeeze(vmask(1,:)').*Vout.west;
            else
              N=size(Uinp.west,3);

              Uout.west=squeeze(0.5.*(Uinp.west(1,:,:)+               ...
                                      Uinp.west(2,:,:)));
              Vout.west=squeeze(0.5.*(Vinp.west(1,1:Mr-1,:)+          ...
                                      Vinp.west(1,2:Mr  ,:)));

              Uout.west=repmat(squeeze(umask(1,:)'),[1,N]).*Uout.west;
              Vout.west=repmat(squeeze(vmask(1,:)'),[1,N]).*Vout.west;
            end
          case {2}
            if (is2d),
              Uout.east=squeeze(0.5.*(Uinp.east(1,:)+                 ...
                                      Uinp.east(2,:)));
              Vout.east=squeeze(0.5.*(Vinp.east(1,1:Mr-1)+            ...
                                      Vinp.east(1,2:Mr  )));

              Uout.east=squeeze(umask(end,:)').*Uout.east;
              Vout.east=squeeze(vmask(end,:)').*Vout.east;
            else
              N=size(Uinp.east,3);

              Uout.east=squeeze(0.5.*(Uinp.east(1,:,:)+               ...
                                      Uinp.east(2,:,:)));
              Vout.east=squeeze(0.5.*(Vinp.east(2,1:Mr-1,:)+          ...
                                      Vinp.east(2,2:Mr  ,:)));

              Uout.east=repmat(squeeze(umask(end,:)'),[1,N]).*Uout.east;
              Vout.east=repmat(squeeze(vmask(end,:)'),[1,N]).*Vout.east;
            end
          case {3}
            if (is2d),
              Uout.south=squeeze(0.5.*(Uinp.south(1:Lr-1,1)+          ...
                                       Uinp.south(2:Lr  ,1)));
              Vout.south=squeeze(0.5.*(Vinp.south(:,1)+               ...
                                       Vinp.south(:,2)));
            
              Uout.south=squeeze(umask(:,1)).*Uout.south;
              Vout.south=squeeze(vmask(:,1)).*Vout.south;
            else
              N=size(Uinp.south,3);

              Uout.south=squeeze(0.5.*(Uinp.south(1:Lr-1,1,:)+        ...
                                       Uinp.south(2:Lr  ,1,:)));
              Vout.south=squeeze(0.5.*(Vinp.south(:,1,:)+             ...
                                       Vinp.south(:,2,:)));
            
              Uout.south=repmat(squeeze(umask(:,1)),[1,N]).*Uout.south;
              Vout.south=repmat(squeeze(vmask(:,1)),[1,N]).*Vout.south;
            end
          case {4}
            if (is2d),
              Uout.north=squeeze(0.5.*(Uinp.north(1:Lr-1,2)+          ...
                                       Uinp.north(2:Lr  ,2)));
              Vout.north=squeeze(0.5.*(Vinp.north(:,1)+               ...
                                       Vinp.north(:,2)));

              Uout.north=squeeze(umask(:,end)).*Uout.north;
              Vout.north=squeeze(vmask(:,end)).*Vout.north;
            else
              N=size(Uinp.north,3);

              Uout.north=squeeze(0.5.*(Uinp.north(1:Lr-1,2,:)+        ...
                                       Uinp.north(2:Lr  ,2,:)));
              Vout.north=squeeze(0.5.*(Vinp.north(:,1,:)+             ...
                                       Vinp.north(:,2,:)));
            
              Uout.north=repmat(squeeze(umask(:,end)),[1,N]).*Uout.north;
              Vout.north=repmat(squeeze(vmask(:,end)),[1,N]).*Vout.north;
            end
        end
      end
    end
  else
    if (is2d),
      Uout=0.5.*(Uinp(1:Lr-1,1:Mr)+Uinp(2:Lr,1:Mr));
      Vout=0.5.*(Vinp(1:Lr,1:Mr-1)+Vinp(1:Lr,2:Mr));

      Uout=umask.*Uout;
      Vout=vmask.*Vout;
    else
      Uout=0.5.*(Uinp(1:Lr-1,1:Mr,:)+Uinp(2:Lr,1:Mr,:));
      Vout=0.5.*(Vinp(1:Lr,1:Mr-1,:)+Vinp(1:Lr,2:Mr,:));

      Uout=repmat(umask,[1,1,N]).*Uout;
      Vout=repmat(vmask,[1,1,N]).*Vout;
    end
  end
end

return
