% create_delilah_2dspec
%
% jcw 22Dec2020
%

% create the offshore 2d spec file in SWAN format from the Delilah data
% at http://www.frf.usace.army.mil/dksrv/del90dir.html
%
cd E:\data\models\InWave\readswan\Projects\Delilah\data\8m_spec
%
%  load the observed spectra from 8m array
%
A=dir;
numfiles=size(A,1);
count=0;
for mm=3:numfiles
  disp(['doing file ',num2str(mm-2),' of ',num2str(numfiles-2),' files.'])
  count=count+1;
  fid=fopen(A(mm).name);
%read time
  line=fscanf(fid,'%f',8);
  yymmdd(count)=floor(line(1)/10000);
  hh(count)=floor(line(1)-yymmdd(count)*10000)/100;
%skip 13 lines
  line=fscanf(fid,'%f',8);
  line=fscanf(fid,'%f',8);
  line=fscanf(fid,'%f',7);
  for nn=1:9
    line=fscanf(fid,'%f',10);
  end
  line=fscanf(fid,'%f',1);
%get freq and scale fac
  for jj=1:29
    line=fscanf(fid,'%f',5);
    freq(jj)=line(2);
    scale=line(3);
    for nn=1:11
      line=fscanf(fid,'%f',8);
      S(count,jj,1+(nn-1)*8:8+(nn-1)*8)=line*scale;
    end
    for nn=12
      line=fscanf(fid,'%f',3);
      S(count,jj,1+(nn-1)*8:3+(nn-1)*8)=line*scale;
    end
  end
  if (line==-1); break; end
  fclose(fid);
end

%
%  ok, now write the data aout.
% 
HEAD=['SWAN   1                                Swan standard spectral file, version'; ...
      '$   Data converted to a SWAN version 40.91A                                 '; ...
      '$   Project: Delilah      ;  run number:                                    '; ...
      'TIME                                    time-dependent data                 '; ...
      '     1                                  time coding option                  '; ...
      'LOCATIONS                               locations in x-y-space              '; ...
      '     1                                  number of locations                 '; ...
      '     950.0000    800.0000                                                   '; ...
      'AFREQ                                   absolute frequencies in Hz          '; ...
      '    29                                  number of frequencies               '];

fid=fopen('8mspec_oct1990.spc2d','w');
for mm=1:size(HEAD,1)
  fprintf(fid,HEAD(mm,:),'%76s \n');
  fprintf(fid,'\n');
end
for mm=1:length(freq)
  fprintf(fid,num2str(freq(mm)),'%s \n');
  fprintf(fid,'\n');
end
fprintf(fid,'NDIR','%s\n');
fprintf(fid,'\n');
fprintf(fid,'    91','%s\n');
fprintf(fid,'\n');
for mm=1:91
  fprintf(fid,num2str(180.0-(mm-1)*2.),'%s \n');
  fprintf(fid,'\n');
end
%
HEAD2=['QUANT                                                                   '; ...
       '     1                                  number of quantities in table   '; ...
       'EnDens                                  energy densities in J/m2/Hz/degr'; ...
       'J/m2/Hz/degr                            unit                            '; ...
       '   -0.9900E+02                          exception value                 '];
for mm=1:size(HEAD2,1)
  fprintf(fid,HEAD2(mm,:),'%72s \n');
  fprintf(fid,'\n');
end
for mm=1:length(yymmdd)
  toadd=['00',num2str(hh(mm))];toadd=toadd(end-1:end);
  datestr=['19',num2str(yymmdd(mm)),'.',toadd,'0000                         date and time'];
  fprintf(fid,datestr,'%s\n');
  fprintf(fid,'\n');
  fprintf(fid,'FACTOR','%s\n');
  fprintf(fid,'\n');
  fprintf(fid,'    1.0E+0','%s\n');
  fprintf(fid,'\n');
  for jj=1:length(freq)
    fprintf(fid,num2str(S(mm,jj,:)*1025.*9.81),'%s \n');
    fprintf(fid,'\n');
  end
end
fclose(fid);


