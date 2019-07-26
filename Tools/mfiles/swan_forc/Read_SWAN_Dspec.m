function [nml,locx,locy,nfreqs,freqhz,ndirs,dir,exc,dateiso,FC,vadens]=Read_SWAN_Dspec(fin);
% function [nml,nfreqs,freqhz,ndirs,dir,exc,dateiso,FC,vadens]=Read_SWAN_Dspec(fin);
% INPUT:
% fin=Input file name
% OUTPUT:
% nml=number of locations
% nfreqs=number of frequencies
% freqhz= frequencies in hertz
% ndirs= number of directions
% dir= Directions in deg
% exc= excepitional value
% dateiso=Date in ISO Format
% FC= Multiplying Factor
% vadens=Spectral Density


fid=fopen(fin);
eofstat = feof(fid);  
for i=1:1:6
    junk=fgetl(fid);
end

C=fgetl(fid);
nml=str2num(C(1:6)); %
for iloc = 1:nml
    C=fgetl(fid);
    locx(iloc)=str2num(C(1:15));   
    locy(iloc)=str2num(C(16:24));  
end

junk=fgetl(fid);    
C=fgetl(fid);
nfreqs=str2num(C(1:6));   %  'number of frequencies'

for i=1:1:nfreqs
C=fgetl(fid);
freqhz(i)=str2num(C(1:10));% 
end


junk=fgetl(fid);     % 'spectral Nautical directions in degr'
C=fgetl(fid);        % 'number of directions'
ndirs=str2num(C(1:6));

for i=1:1:ndirs
    C=fgetl(fid); 
    dir(i)=str2num(C(1:10)); 
end

j=find(dir<0);
dir(j)=dir(j)+360;

junk=fgetl(fid);      %    WRITE (9, 102) 'QUANT'
C=fgetl(fid); 
numq=str2num(C(1:6)); % 'number of quantities in table'
junk=fgetl(fid);      % 'VaDens', 'variance densities in m2/Hz/degr'
junk=fgetl(fid);      % 'm2/Hz/degr',  'unit'
C=fgetl(fid); 
exc=str2num(C(1:6));  % 'exception value'

index=0;

while (feof(fid)==0)
   index=index+1;
   C=fgetl(fid); 
   dateiso(index)=str2num(C(1:16)); % 'Time in ISO format' 
   disp(C(1:16))
   for iloc =1:nml
       junk=fgetl(fid);                % 'FACTOR', 'multiplication factor'
       C=fgetl(fid);
       FC(index,iloc)=str2num(C(1:18));
       for ifr = 1:nfreqs
           C=fgetl(fid);
           vadens(ifr,1:ndirs,index,iloc)=strread(C,' %10.5f',ndirs);
       end
   end
end
