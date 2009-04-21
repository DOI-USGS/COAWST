function [s]=swan_meta_v3(ifile);
% SWAN_META returns swan meta information
% Usage: [s]=swan_meta_v2(ifile);
%
% Input:   ifile = swan input file name (e.g. 'INPUT');
%
% Output:  s.xsize           = x grid size
%          s.ysize           = y grid size
%          s.jdmat0          = time of first output field (DATENUM format) 
%          s.jdmat1          = time of end output field (DATENUM format) 
%          s.dt              = DT between output (seconds)
% These added by Brandy Armstrong 6/16/08 for xml docs on THREDDS 
%          s.gen             = Komen,Janssen or Westhuysen
%          s.bins            = number of bins 
%          s.var             = variable names
%          s.matfile         = matfile names
%
% Example: s=swan_meta('INPUT');

% Note: It is assumed that strings 'CGRID' and 'UBOT' will be found, as these
% are used to locate the grid size and the output time information

% Rich Signell (rsignell@usgs.gov)
% 26-Sep-2005
% 27-Sep-2005 jcwarner, added catch for multiplier index after UBOT.
% hsk 06-Oct-2006: Added to pick (start,end,intervals) of simulated output
%                   records. 
%                       start-time: s.jdmat0
%                       end-time:   s.jdmat1
%                       interval:   s.dt
% hsk 26-Dec-2006: added to pick # of records (nt) in each simulation
%                  output data
%                       # of records:  s.nt

rehash

f=textread(ifile,'%s');
%fid = fopen(ifile);

% find xsize,ysize
ind=find(strcmp(f,'CGRID'));
if isempty(ind),
  disp('Error: could not determine size of grid')
  return
else
  temp=char(f(ind+1));
  if temp(1:7)=='REGULAR'
    s.xsize=str2num(char(f(ind+7)))+1;
    s.ysize=str2num(char(f(ind+8)))+1;
  else
    s.xsize=str2num(char(f(ind+2)))+1;
    s.ysize=str2num(char(f(ind+3)))+1;
  end
end

% find start and end times 
indx=find(strcmp(f,'COMPUTE'));
tmp=char(f(indx+2));        % lookfor 
tmp2=char(f(indx+1));
%if (length(tmp)==15 & length(char(f(indx+5)))==15)
if (length(tmp)==15 & strcmp(tmp2(1:3),'NON'))
    g1=tmp;                 % start-time
    g2=char(f(indx+5));     % end-time
elseif (length(tmp)==15 & strcmp(tmp2(1:3),'STA'))
    g1=tmp;                 % start-time
    g2=g1;                  % end-time
    disp('set end time = to start time')
else
    error('Error: could not determine OUTPUT start time')
    return
end

G1=iso8601(g1);
G2=iso8601(g2);

s.jdmat0=G1(1);
s.jdmat1=G2(1);

% find DT
ind=find(strcmp(f,'HSIGN'),1); %had to add 1 here to find only first occurance BNA
tmp=[char(f(ind+1)), '   ']; %jcw
%If next string is the scalar multiplier then add 1 to index. jcw
if ~strcmp(upper(tmp(1:3)),'OUT')
  ind=ind+1;
  tmp=[char(f(ind+1)), '   '];
end
if strmatch(upper(tmp(1:3)),'OUT')  % check to make sure this is an OUTPUT line
  tmp=char(f(ind+2));
  if length(tmp)==15,  % check if length is correct for an ISO time string, if not, take next
%    g=tmp;
    dt=str2num(char(f(ind+3)));
    dtunits=char(f(ind+4));
  elseif length(char(f(ind+3)))==15, 
%    g=char(f(ind+3));
    dt=str2num(char(f(ind+4)));
    dtunits=char(f(ind+5));
  else
    disp('Error: could not determine OUTPUT time interval')
    return
  end
  switch dtunits
    case 'DAY'
    s.dt=dt*24*3600; 
    case 'HR'
    s.dt=dt*3600; 
    case 'MIN'
    s.dt=dt*60;
    case 'SEC'
    s.dt=dt;
    otherwise
    disp('Error: could not determine OUTPUT start time and time interval')
    return
  end
else
  disp('Error: could not determine OUTPUT start time and time interval')
end

% find GEN swan.gen
ind=find(strcmp(f,'GEN3'));
if isempty(ind),
    disp('Error: could not determine GEN3')
    return
else
    temp=char(f(ind+1));
    if strcmp(temp(1:4),'PROP')
        s.gen='KOMEN';
    elseif sctrcmp(temp(1:5),'KOMEN')
        s.gen='KOMEN';
    elseif strcmp(temp(1:5),'WESTH')
        s.gen='WESTHUYSEN';
    elseif sctrcmp(temp(1:7),'JANSSEN')
        s.gen='JANSSEN';
    end
end

% find bins swan.bins
ind=find(strcmp(f,'CIRCLE'));
if isempty(ind),
    disp('Error: could not determine BINS')
    return
else
    temp=char(f(ind+1));
    s.bins=str2num(temp);
end

% find name of run swan.name
ind=find(strcmp(f,'PROJECT'));
if isempty(ind),
    disp('Error: could not determine PROJECT')
    return
else
    temp=char(f(ind+1));
    temp2=char(f(ind+2));
    s.name=[temp temp2];
end


%find variables swan.var and output s.matfile
ind=find(strcmp(f,'NOHEADER'));
[ir ic]=size(ind);
for i=1:ir
    temp=char(f(ind(i)+1));
    s.var{i}=(temp(2:end-5));
    s.matfile{i}=(temp);
end

%find waves input
ind=find(strcmp(f,'BOUNDSPEC'));
[ir ic]=size(ind);
    temp=char(f(ind+10));
    spec=strcmp(temp(1,2:end-2),['../' s.name(end-5:end-1) '/POINT1.spc2d']);
    tpar=strcmp(temp(1,2:end-2),['../forcing/TPAR1.txt']);
    if tpar
        s.spec='Wave Watch 3 TPAR';
    elseif spec
        s.spec='SWAN 2d spectra';
    end
    s.specnum=ir;





s.wind='Rotated combined and smoothed NARR and HRD data';
% s.spec='Wave Watch 3 TPAR files';
s.input=ifile;
% find NT (hsk 12/26/06)
nd=(s.jdmat1-s.jdmat0);      % # of days for records
s.nt=(nd*24*3600)/s.dt + 1;

f2=textread('run_swan','%s');
ind=find(strcmp(f2,'cd'));
temp=char(f2(ind+1));
s.id=temp(end-17:end);

save swan


%fclose(fid);
