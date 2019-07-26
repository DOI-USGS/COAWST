function createswan2Dspec(YYYYMODD,HHMMSS,Loc,Freq,Dir,Spec,Fname)
%function createswan2Dspec(YYYYMODD,HHMMSS,Loc,Freq,Dir,Spec)
%USAGE: This function takes time series of 2D spectrum from wavewatch/instrument to create a two-dimensional
%spectral input file. 

%INPUT: 
%YYYYMODD= Part of time string corresponding to year, month and day (e.g., 20100201)
%HHMMSS  = Part of time string corresponding to hour, min and ss (e.g., 064500)
%Loc     = Co-ordinates for location at which boundary information is provided (i.e., Loc = [Lon,Lat])
%Freq    = Vector of frequency spacing of the two-dimensional spectrum
%Dir     = Vector of directional spacing of the two-dimensional spectrum
%Spec    = Frequency-directional spectrum (expected units is m^2/Hz/Deg). Also it is 
%          assumed that the spectrum is a 4-dimensional variable with a structure 
%          Spec(NT,NLOC,NFREQ,NDIR). 
%Fname   = File name for output 

if (nargin<7 || isempty(Fname)==1)
	Fname= 'wave_forcing.bnd'
end 

if nargin<5
	error('Function requires time index, location, frequency, direction and 2-D spec')
end

NT    = length(YYYYMODD);     %Number of time steps
NLOC  = size(Loc,1)     ;     %Number of locations
NFREQ = length(Freq)    ;     %Number of frequencies
NDIR  = length(Dir)     ;     %Number of directions
rho   = 1025            ;     %SWAN allows changing rho value, but deafult is 1025 kg/m^3
g     = 9.81            ;     %Accn. due to gravity
ver   = '40.85'			;     %SWAN version
proj  = 'cape_run'      ;	  %Project Name
tcode = 1				;     %Time coding option
Exval = -99				;	  %Exception value
Time  = YYYYMODD+HHMMSS/(1000000);  %For correct time format in file

fid   = fopen(Fname,'w');
fprintf(fid,'%s\n','SWAN   1                                Swan standard spectral file, version');
fprintf(fid,'%s\n',['$   Data produced by SWAN version ',ver]);
fprintf(fid,'%s\n',['$   Project: ',proj,'        ;  run number:']);
fprintf(fid,'%s\n','TIME                                    time-dependent data');
fprintf(fid,'%6i',tcode);
fprintf(fid,'%s\n','                                  time coding option');
fprintf(fid,'%s\n','LONLAT                                  locations in spherical coordinates');
fprintf(fid,'%6i',NLOC);
fprintf(fid,'%s\n','                                  number of locations');
for i=1:1:NLOC
	fprintf(fid,'%12.6f %12.6f\n',Loc(i,:));
end
clear i

fprintf(fid,'%s\n','AFREQ                                   absolute frequencies in Hz');
fprintf(fid,'%6i',NFREQ);
fprintf(fid,'%s\n','                                  number of frequencies');
for i=1:1:NFREQ
	fprintf(fid,'%10.4f\n',Freq(i));
end
clear i
fprintf(fid,'%s\n','NDIR                                    spectral nautical directions in degr);');
fprintf(fid,'%6i',NDIR);
fprintf(fid,'%s\n','                                  number of directions');
for i=1:1:NDIR
    if Dir(i)>360
        Dir(i)=Dir(i)-360;
    end
	fprintf(fid,'%10.4f\n',Dir(i));
end
clear i
fprintf(fid,'%s\n', 'QUANT');
fprintf(fid,'%s\n', '     1                                  number of quantities in table');
fprintf(fid,'%s\n', 'VaDens                                  variance densities in m2/Hz/degr');
fprintf(fid,'%s\n', 'm2/Hz/degr                             unit');
%fprintf(fid,'%s\n', 'EnDens                                   variance densities in m2/Hz/degr');
%fprintf(fid,'%s\n', 'J/m2/Hz/degr                             unit');
fprintf(fid,'%14.4e',Exval);
fprintf(fid,'%s\n', '                          exception value');

for p=1:1:NT
    fprintf(fid,'%15.6f',Time(p));
	fprintf(fid,'%s\n','                         date and time');
	for q=1:1:NLOC
        %Edens = rho*g*squeeze(Spec(p,q,:,:));
        Edens = squeeze(Spec(p,q,:,:));
        [ID]=max(max(Edens));
        if ID>(10^-10)
            FAC=(1.01*ID*10^(-4));
        else
            FAC=10^(-5);
        end
        Edens = Edens/FAC;
        FACS  = find_exp(FAC);
		fprintf(fid,'%s\n','FACTOR');
% 		fprintf(fid,'%18.8E\n',FAC);
		fprintf(fid,'%s\n',FACS);
		for r=1:1:NFREQ
            fprintf(fid,'%6d',round(Edens(r,:)));
			fprintf(fid,'\n');
		end
	end
end
			
fclose(fid);
end

function [N]=find_exp(x);
% 
% Convert number x into a string with an exponent format
% xx.xxxxxxxxE-nn

i=-1;
if floor(x)==0
    g=0;
    while g<1
        i=i+1;
        g=floor(rem(x*(10^i),10));
    end
    N = -i;
    G=x*(10^-N);
else
    g=500;
    while (g~=0)
        i=i+1;
        g=floor(x/(10^i));
    end
    N = i-1;
    G=x*(10^-N);
end
N=sprintf('%13.8f%s%+03i',G,'E',N);
end



