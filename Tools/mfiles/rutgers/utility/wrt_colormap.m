function wrt_colormap(Cname, Ncolors)

% WRT_COLORMAP:  Writes colormap into an ASCII file
%
% wrt_colormap(Cname, Ncolors)
%
% This function writes colormap RGB triplets into an ASCII file,
% which it is used to generate color palettes for ROMS plotting
% package.
%
% On Input:
%
%    Cname         Colormap name (string).
%
%    Ncolors       Number of colors to process (scalar or vector)
%
% Example:
%
%    wrt_colormap('viridis', [60 90 120 255])
%
% or
%
%    wrt_colormap('-viridis', 30)
%
% for inverse palette, flipud(Cname).
%
% It uses the "cmocean.m" function to extract the colormap for the
% specified number of colors.
%
% Then, the program GenPal.F can be use to generate ROMS plotting
% package color palette.
%

% svn $Id: wrt_colormap.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
  
%--------------------------------------------------------------------------
% Write out colormap.
%--------------------------------------------------------------------------

for i=1:length(Ncolors),

% Extract colormap.

  Cmap=cmocean(Cname, Ncolors(i));
  [Im,Jm]=size(Cmap);
  
% Write out colormap.

  if (Cname(1:1) == '-'), 
    Fname=strcat(Cname(2:end), '_', num2str(Ncolors(i)), '_flip.dat');
  else
    Fname=strcat(Cname, '_', num2str(Ncolors(i)), '.dat');
  end

  fid=fopen(Fname,'w');

  for i=1:Im,
    R=Cmap(i,1);
    G=Cmap(i,2);
    B=Cmap(i,3);
    fprintf(fid, '%23.16e  %23.16e %23.16e\n', R, G, B);
  end

  fclose(fid);
end

return
