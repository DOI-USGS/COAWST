function [freq,dir_sort,E2]=read_swan_spc(filename)

fid=fopen(filename,'r');
for nn=1:10
    tline = fgetl(fid)
end
nfreq=str2num(tline(5:6));
freq=zeros(nfreq,1);
for ff=1:nfreq
    tline = fgetl(fid);
    freq(ff)=str2num(tline(1:end));
end
tline = fgetl(fid);
tline = fgetl(fid);
ndir=str2num(tline(5:6));
dir_sort=zeros(ndir,1);
for dd=1:ndir
    tline = fgetl(fid);
    dir_sort(dd)=str2num(tline(1:end));
end

for ii=1:8
    tline = fgetl(fid);
end
Sfact=str2num(tline(5:end));

E2 = fscanf(fid,'%i',[ndir,nfreq])';
E2=E2.*Sfact;
E2=E2/1025/9.81;   % convert Joules/m2 to m2

fclose(fid)

end

