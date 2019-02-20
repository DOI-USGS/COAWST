function [umask,vmask,pmask]=uvp_masks(rmask)

%
% UVP_MASK:  Computes Land/Sea U-, V- and PSI-mask
%
% [umask,vmask,pmask]=uvp_masks(rmask)
%
%  This function computes the Land/Sea masks on U-, V-, and PSI-points
%  from the mask on RHO-points.
%
%  On Input:
%
%    rmask        Land/Sea mask on RHO-points (real matrix).
%
%  On Output:
%
%    umask        Land/Sea mask on U-points (real matrix).
%    vmask        Land/Sea mask on V-points (real matrix).
%    pmask        Land/Sea mask on PSI-points (real matrix).
%

% svn $Id: uvp_masks.m 916 2018-07-14 01:28:47Z arango $
%=========================================================================%
%  Copyright (c) 2002-2018 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arang       %
%=========================================================================%

[Lp,Mp]=size(rmask);
L=Lp-1;
M=Mp-1;

%  Land/Sea mask on U-points.

umask(1:L,1:Mp)=rmask(2:Lp,1:Mp).*rmask(1:L,1:Mp);

%  Land/Sea mask on V-points.

vmask(1:Lp,1:M)=rmask(1:Lp,2:Mp).*rmask(1:Lp,1:M);

%  Land/Sea mask on PSI-points. The old computation of PSI-mask is replaced
%  with the newer scheme.  Notice the shift because Matlab doesnot support
%  zero index arrays.
%
%  pmask(1:L,1:M)=rmask(1:L,1:M ).*rmask(2:Lp,1:M ).*                   ...
%                 rmask(1:L,2:Mp).*rmask(2:Lp,2:Mp);

pmask=nan([L M]);

for jr=2:Mp
  jp=jr-1;
  for ir=2:Lp
    ip=ir-1;
    if ((rmask(ir-1,jr  ) > 0.5) &&                                     ...
        (rmask(ir  ,jr  ) > 0.5) &&                                     ...
        (rmask(ir-1,jr-1) > 0.5) &&                                     ...
        (rmask(ir  ,jr-1) > 0.5))
      pmask(ip,jp)=1;
    elseif ((rmask(ir-1,jr  ) < 0.5) &&                                 ...
            (rmask(ir  ,jr  ) > 0.5) &&                                 ...
            (rmask(ir-1,jr-1) > 0.5) &&                                 ...
            (rmask(ir  ,jr-1) > 0.5))
      pmask(ip,jp)=1;
    elseif ((rmask(ir-1,jr  ) > 0.5) &&                                 ...
            (rmask(ir  ,jr  ) < 0.5) &&                                 ...
            (rmask(ir-1,jr-1) > 0.5) &&                                 ...
            (rmask(ir  ,jr-1) > 0.5))
      pmask(ip,jp)=1;
    elseif ((rmask(ir-1,jr  ) > 0.5) &&                                 ...
            (rmask(ir  ,jr  ) > 0.5) &&                                 ...
            (rmask(ir-1,jr-1) < 0.5) &&                                 ...
            (rmask(ir  ,jr-1) > 0.5))
      pmask(ip,jp)=1;
    elseif ((rmask(ir-1,jr  ) > 0.5) &&                                 ...
            (rmask(ir  ,jr  ) > 0.5) &&                                 ...
            (rmask(ir-1,jr-1) > 0.5) &&                                 ...
            (rmask(ir  ,jr-1) < 0.5))
      pmask(ip,jp)=1;
    elseif ((rmask(ir-1,jr  ) > 0.5) &&                                 ...
            (rmask(ir  ,jr  ) < 0.5) &&                                 ...
            (rmask(ir-1,jr-1) > 0.5) &&                                 ...
            (rmask(ir  ,jr-1) < 0.5))
      pmask(ip,jp)=2;
    elseif ((rmask(ir-1,jr  ) < 0.5) &&                                 ...
            (rmask(ir  ,jr  ) > 0.5) &&                                 ...
            (rmask(ir-1,jr-1) < 0.5) &&                                 ...
            (rmask(ir  ,jr-1) > 0.5))
      pmask(ip,jp)=2;
    elseif ((rmask(ir-1,jr  ) > 0.5) &&                                 ...
            (rmask(ir  ,jr  ) > 0.5) &&                                 ...
            (rmask(ir-1,jr-1) < 0.5) &&                                 ...
            (rmask(ir  ,jr-1) < 0.5))
      pmask(ip,jp)=2;
    elseif ((rmask(ir-1,jr  ) < 0.5) &&                                 ...
            (rmask(ir  ,jr  ) < 0.5) &&                                 ...
            (rmask(ir-1,jr-1) > 0.5) &&                                 ...
            (rmask(ir  ,jr-1) > 0.5))
      pmask(ip,jp)=2;
    else
      pmask(ip,jp)=0;
    end
  end
end

return
