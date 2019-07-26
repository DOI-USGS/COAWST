function [var_rho]=v2rho_2d_mw(var_v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% transfert a field at u points to a field at rho points
%
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2001-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Lp,M]=size(var_v);
Mp=M+1;
Mm=M-1;
var_rho=zeros(Lp,Mp);
var_rho(:,2:M)=0.5*(var_v(:,1:Mm)+var_v(:,2:M));
var_rho(:,1)=var_rho(:,2);
var_rho(:,Mp)=var_rho(:,M);

return

