function [FA,DIR,S1]=read_swan(filename) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%          READS SWAN 2D SPECTRA     %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% open the file for reading
fid = fopen(filename,'r');

done=0;
%scan file to get number of freqs
while (~done)
  tline=fgetl(fid);
  if strncmp(tline,'AFREQ',5)
    tline=fgetl(fid);
    numfreq=str2num(tline(1:20));
    fa = fscanf(fid, '%g', [1 numfreq]);
  end
  if strncmp(tline,'NDIR',4)
    tline=fgetl(fid);
    numdir=str2num(tline(1:20));
    dir = fscanf(fid, '%g', [1 numdir]);
  end
  if strncmp(tline,'FACTOR',6)
    fact = fscanf(fid, '%g', [1]);
    S = fscanf(fid, '%g', [ numdir numfreq]);
    S = fact.*S.';
    done=1;
  end
end
fclose(fid)


[FA,DIR]=meshgrid(fa,dir);
FA=FA';
DIR=DIR';
S1=S;

