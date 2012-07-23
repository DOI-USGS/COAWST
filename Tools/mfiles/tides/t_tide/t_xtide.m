function pred=t_xtide(varargin);
% T_XTIDE Tidal prediction
% YOUT=T_XTIDE(STATION) makes a tidal prediction
% for the current day using the harmonics file from XTIDE. 
% if STATION is a string then the first match found in the database is
% used, you can request matches to other stations by appending '(2)' to
% the string. If you don't know the station name but want to find the
% nearest to a given LONG,LAT, try T_XTIDE(LONG,LAT).
%
% The times of predicted tides are given by the next numerical argument
% (if any), e.g. [...]=T_XTIDE(STATION,TIM). 
% TIM can be: a vector of matlab-format decimal days (from DATENUM).
%           : a scalar <1000, taken as the number of days from present
%           : a scalar >1000, taken as the starting time in matlab-format 
%             for a 2 day time series. 
%           : not given, in which case the current time is used as a start 
%             time.
%
% Times are usually taken to be in 
% standard time for the given location (no daylight savings adjustments); 
% if in doubt use the 'info' or 'full' options where offset from UTC is given.
%
% 
% Other optional arguments can be specified following this using 
% property/value pairs: 
%
%     'format'     'raw' (default)
%                    YOUT is just a time series to match the time in TIM
%
%                  'times'
%                    YOUT is a structure of high/low tide information
%                    between times min(TIM) and max(TIM).
%
%                  'info'
%                    YOUT is a structure giving station information
%                    (location, time zone, units, harmonic constituents)
%
%                  'full'
%                    Combination of 'raw' and 'info' in a structure YOUT.
%
%     'units'     {'meters' | 'feet' | 'm/s' | 'knots' | 'original' }
%                    Units of result (default is original units)
%
% If no output argument is specified data is plotted and/or displayed.
%
% If the chosen name matches several stations then the first in the list is
% chosen. Specific choices can be made by appending a '(2)' or '(3)' (etc.)
% to the name, e.g.  T_XTIDE('tofino (2)',...).
%
%  Requires the xtide harmonics file  - get this from 
%            http://bel-marduk.unh.edu/xtide/files.html

% R. Pawlowicz 1/Dec/2001
% Version 1.0
%          16/May/02 - added lat/long options (thanks to Richard Dewey).


% Get the harmonics data from a) a mat-file if it exists, or b) from a harmonics
% file.

if ~exist('t_xtide.mat','file'), % Read the harmonics file and make a mat file

  filnam='/usr/share/xtide/harmonics.txt';
  
  fprintf('\n********Can''t find mat-file t_xtide.mat ********\n\n');
  fprintf('Attempting to generate one from an xtide harmonics file....\n\n');
  fprintf('Latest version available from http://bel-marduk.unh.edu/xtide/files.html\n\n');
  
  % Input name
  fid=-1;
  while fid==-1,
    rep=filnam;
    while (lower(rep(1))~='y'),
     filnam=rep;
     rep='n';
     rep=input(['Harmonics filename: ' filnam '? (y/Y/new file name):'],'s');
     if isempty(rep), rep='y'; end;
    end; 
    
    fid=fopen(filnam);
    if fid==-1,
      fprintf(['\n****** Can''t open filename ->' filnam '<-\n\n']);
    end;
  end;
    
  fprintf('Reading harmonics file (this will take a while)\n');
  [xtide,xharm]=read_xtidefile(fid);
  
  fprintf('Saving harmonic information to t_xtide.mat\n');
  save t_xtide xtide xharm
   
else
  load t_xtide
end;

if nargin>0,

  if isstr(varargin{1}),  % Station  name given
    % Identify station - look for exact match first
    ista=strmatch(lower(varargin{1}),lower(xharm.station),'exact');
    % otherwise go for partial matches
    if isempty(ista),
      % First check to see if a number was selected:
      inum=-10;
      while inum<-1,
        inum=inum+1;
        ll=findstr(lower(varargin{1}),sprintf('(%d)',-inum));
        if ~isempty(ll),
          inum=abs(inum);
	      varargin{1}=deblank(varargin{1}(1:ll-1));
        end;
      end;  
      ista=strmatch(lower(varargin{1}),lower(xharm.station));
      if length(ista)>1,
        if inum>0 & inum<=length(ista),
          ista=ista(inum);
        else	
          fprintf('Ambiguous Station Choice - Taking first of:\n');
          for kk=1:length(ista),
	        fprintf('%5d: %s\n',ista(kk),deblank(xharm.station(ista(kk),:)));
	        fprintf('      Long: %.4f  Lat: %.4f \n',xharm.longitude(ista(kk)),xharm.latitude(ista(kk)));
          end;
          fprintf('\n');
          ista=ista(1);
        end 	
      elseif length(ista)==1 & inum>1,
        fprintf('***Can''t find variant (%d) of station - Taking only choice\n',inum);
      elseif length(ista)==0,  
        error('Could not match station');
      end;    
     end;
     varargin(1)=[];

   else   % Lat/long?
      [dist,hdg]=t_gcdist(xharm.latitude,xharm.longitude,varargin{2},varargin{1});
      [mind,ista]=min(dist);
      if length(ista)>1,
        fprintf('Ambiguous Station Choice - Taking first of:\n');
        for kk=1:length(ista),
	      fprintf('%5d: %s\n',ista(kk),deblank(xharm.station(ista(kk),:)));
	      fprintf('      Long: %.4f  Lat: %.4f \n',xharm.longitude(ista(kk)),xharm.latitude(ista(kk)));
        end;
        fprintf('\n');
        ista=ista(1);
      else
 	    fprintf('%5d: %s\n',ista,deblank(xharm.station(ista,:)));
	    fprintf('      Long: %.4f  Lat: %.4f \n',xharm.longitude(ista),xharm.latitude(ista)); 
      end;
      varargin(1:2)=[];
   end;
  
  % Time vector (if available) otherwise take current time.

  if length(varargin)>0 & ~isstr(varargin{1}),
    tim=varargin{1};
    tim=tim(:)';
    varargin(1)=[];
    if length(tim)==1,
      if tim<1000,
        dat=clock;
        tim=datenum(dat(1),dat(2),dat(3))+[0:1/48:tim];
      else
        tim=tim+[0:1/48:2]; % 2 days worth.
      end;	 	
    end;
  else 
    dat=clock;
    tim=datenum(dat(1),dat(2),dat(3))+[0:.25:48]/24;
  end;
 
   % Parse properties

  format='raw';
  unt='original';
  
  k=1;
  while length(varargin)>0,
      switch lower(varargin{1}(1:3)),
	case 'for',
	 format=lower(varargin{2});
	case 'uni',
	 unt=lower(varargin{2}); 
	otherwise,
           error(['Can''t understand property:' varargin{1}]);
      end;
      varargin([1 2])=[]; 
  end;
 
  % if we want a time series
  pred=[];
  % Convert units if requested.
  [units,convf]=convert_units(unt,xharm.units(ista,:));
  if strcmp(format(1:2),'ra') | strcmp(format(1:2),'fu') | strcmp(format(1:2),'ti')
    
    % Data every minute for hi/lo forecasting.
    if strcmp(format(1:2),'ti'),
      tim=tim(1):(1/1440):tim(end); 
    end;

    % Convert into time since the beginning of year
    mid=datevec(mean(tim));
    iyr=mid(1)-xtide.startyear+1;
    lt=length(tim);
    xtim=(tim-datenum(mid(1),1,1))*24; % Hours since beginning of year

    %-----------------------------------------------------
    % Sum up everything for the prediction!

    pred=xharm.datum(ista)+sum( ...
      repmat(xtide.nodefactor(:,iyr).*xharm.A(ista,:)',1,lt).* ...
      cos( ( xtide.speed*xtim + repmat(xtide.equilibarg(:,iyr)-xharm.kappa(ista,:)',1,lt) )*(pi/180) ),1);
    %-----------------------------------------------------

    pred=pred*convf;
    
    % Compute times of hi/lo from every-minute data
    if strcmp(format(1:2),'ti'),
     % Check if this is a current station
      if ~isempty(findstr('Current',xharm.station(ista,:))), currents=1; else currents=0; end;
      dpred=diff(pred);
      ddpred=diff(dpred>0);

      flat=find(ddpred~=0)+1;
      slk=find(sign(pred(1:end-1))~=sign(pred(2:end)));
      
      hi.mtime=tim(flat);
      hi.value=pred(flat);

      hi.type=zeros(size(flat));
      hi.type(find(ddpred(flat-1)<0))=1;  % 0=lo, 1=hi
      hi.units=deblank(units);
      
      pred=hi;
    end;
  end;
  
  % Create information structure
  
  if strcmp(format(1:2),'in') | strcmp(format(1:2),'fu'),
    if ~isempty(pred), 
      pred.yout=pred; 
      pred.mtime=tim; 
    else
      kk=find(xharm.A(ista,:)~=0);
      pred.freq=xtide.name(kk,:);
      pred.A=full(xharm.A(ista,kk)')*convf;
      pred.kappa=full(xharm.kappa(ista,kk)'); 
    end;
    pred.station=deblank(xharm.station(ista,:));
    pred.longitude=xharm.longitude(ista);
    pred.latitude=xharm.latitude(ista);
    pred.timezone=xharm.timezone(ista);
    pred.units=deblank(units);
    pred.datum=xharm.datum(ista)*convf;
  end;
 
end;

% If no output parameters then we plot or display things

if nargout==0,
  switch format(1:2),    
    case 'ti',
  
    fprintf('High/Low Predictions for %s\n',xharm.station(ista,:));
    fprintf('Time offset %.1f from UTC\n\n',xharm.timezone(ista));
    
    outstr=repmat(' ',length(flat),41);
    outstr(:,1:20)=datestr(hi.mtime);
    outstr(:,22:27)=reshape(sprintf('%6.2f',hi.value),6,length(flat))';
    if currents,
      ll=hi.type==1;
      outstr(ll,31:41)=repmat(' Flood Tide',sum(ll),1);
      ll=hi.type==0;
      outstr(ll,31:41)=repmat(' Ebb Tide  ',sum(ll),1);
    else
      ll=hi.type==1;
      outstr(ll,31:41)=repmat(' High Tide ',sum(ll),1);
      ll=hi.type==0;
      outstr(ll,31:41)=repmat(' Low Tide  ',sum(ll),1);
    end;
    disp(outstr)   
           
    case 'ra'
     plot(tim,pred)
     datetick;
     title(['Tidal prediction for ',deblank(xharm.station(ista,:)) ' beginning ' datestr(tim(1))]); 
     ylabel(deblank(xharm.units(ista,:)));

    case 'fu'
     plot(tim,pred.yout)
     datetick;
     title(['Tidal prediction for ',deblank(xharm.station(ista,:)) ' beginning ' datestr(tim(1))]); 
     ylabel(deblank(xharm.units(ista,:)));
     
    case 'in',
        
    fprintf('Station: %s\n',pred.station);
    if pred.longitude<0, lon='W'; else lon='E'; end;
    if pred.latitude<0,  lat='S'; else lat='N'; end;
    fprintf('Location: %d %.1f'' %c, %d %.1f'' %c\n',fix(abs(pred.latitude)),rem(abs(pred.latitude),1)*60,...
         lat,fix(abs(pred.longitude)),rem(abs(pred.longitude),1)*60,lon);
    fprintf('Time offset %.1f from UTC\n\n',pred.timezone);
     	
   end;
   clear pred
end;  
  
%%%%%%%%%%%%%%%%%%%%
function [xtide,xharm]=read_xtidefile(fid);
% Reads the xtide harmonics file and creates a data structure
% with all that info for faster access


l=fgetl_nocom(fid);

ncon=sscanf(l,'%d');

xtide=struct('name',repmat(' ',ncon,8),'speed',zeros(ncon,1),...
	     'startyear',0,'equilibarg',zeros(ncon,68),'nodefactor',zeros(ncon,68));

for k=1:ncon,
 l=fgetl_nocom(fid);
 xtide.name(k,:)=l(1:8);
 xtide.speed(k)=sscanf(l(9:end),'%f');
end;

xtide.startyear=sscanf(fgetl_nocom(fid),'%d');

nyear=sscanf(fgetl_nocom(fid),'%d');

for k=1:ncon,
  l=fgetl(fid);
  xtide.equilibarg(k,:)=fscanf(fid,'%f',nyear);
  l=fgetl(fid);
end;
l=fgetl(fid); % Skip *END*

nyear=sscanf(fgetl_nocom(fid),'%d');

for k=1:ncon,
  l=fgetl(fid);
  xtide.nodefactor(k,:)=fscanf(fid,'%f',nyear);
  l=fgetl(fid);
end;
l=fgetl(fid); % Skip *END*

% Now read in all harmonic data


%nsta=1754; % This is number of stations in harmonics (1998-07-18)
%nsta=3351; % This is number of stations in v1.42 or harmonics file
nsta=3316; % This is number in v1.51

xharm=struct('station',repmat(' ',nsta,79),'units',repmat(' ',nsta,8),...
	     'longitude',zeros(nsta,1),'latitude',zeros(nsta,1),...
	     'timezone',zeros(nsta,1),'datum',zeros(nsta,1),...
	     'A',zeros(nsta,ncon),'kappa',zeros(nsta,ncon));

nh=0;
while length(l)>0 & l(1)~=-1,
 
  l=[l '   '];
  nh=nh+1;
  while ~strcmp(l(1:3),'# !'),
    l=[fgetl(fid) '   '];
  end;
  while strcmp(l(1:3),'# !'),
   switch l(4:7),
    case 'unit',
     tmp=deblank(l(findstr(l,':')+2:end));
     xharm.units(nh,1:length(tmp))=tmp;
    case 'long',
      xharm.longitude(nh)=sscanf(l(findstr(l,':')+1:end),'%f');
    case 'lati'  
      xharm.latitude(nh)=sscanf(l(findstr(l,':')+1:end),'%f');
   end;
   l=fgetl(fid);
  end; 
  tmp=deblank(l);
  if tmp(1)~='#', % Not commented out
    xharm.station(nh,1:length(tmp))=tmp;

    tmp=fgetl(fid);
    k=min(findstr(tmp,':'));
    tim=sscanf(tmp(1:k-1),'%d')+sscanf(tmp(k+[1:2]),'%d')/60;
    xharm.timezone(nh)=tim;

    tmp=fgetl(fid);
    xharm.datum(nh)=sscanf(tmp,'%f');

    for k=1:ncon,
      l=fgetl(fid);
      if l(1)~='x',
	ll=min([findstr(' ',l) find(abs(l)==9)]); % space or tab
	tmp=sscanf(l(ll+1:end),'%f',2);
	xharm.A(nh,k)=tmp(1);
	xharm.kappa(nh,k)=tmp(2);
      end;
    end;
    l=fgetl(fid);
  else
    nh=nh-1;  
  end;
  
  if rem(nh,50)==0, fprintf('.'); end;
end;
fprintf('\n');

% Convert internally to sparse matrix storage (much smaller).
xharm.A=sparse(xharm.A);
xharm.kappa=sparse(xharm.kappa);

return;
  
%%%%%%%%%%%%%%%%%%%%  
function l=fgetl_nocom(fid);
% Gets a line that isn't a comment line
%
l=fgetl(fid);
while length(l)>0 & l(1)=='#',
  l=fgetl(fid);
end;
  
%%%%%%%%%%%%%%%%%%%%%
function [units,convf]=convert_units(unt,origunits);
% Conversion factors from origianl units if requested and possible
% (no conversions from knots to feet).
%
  if strcmp(unt(1:3),origunits(1:3)) | strcmp(unt(1:3),'ori'),
    units=origunits;
    convf=1;
  else
   switch unt(1:3),
    case 'fee',
       if strcmp(origunits(1:3), 'met'),
  	  units='feet';
  	  convf=3.2808399;
       else
  	  units=origunits;
  	  convf=1;
       end;
    case 'met',
       if strcmp(origunits(1:3), 'fee'),
  	  units='meters';
  	  convf=0.3048;
       else
  	  units=origunits;
  	  convf=1;
       end;
    case 'm/s',
       if strcmp(origunits(1:3), 'kno'),
  	  units='meters/sec';
  	  convf=0.51444444;
       else
  	  units=origunits;
  	  convf=1;
       end;
    case 'kno',
       if strcmp(origunits(1:3), 'm/s'),
  	  units='knots';
  	  convf=1.9438445;
       else
  	  units=origunits;
  	  convf=1;
       end;
    otherwise
      error('Unknown units')
    end;
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d,hdg]=t_gcdist(lat1,lon1,lat2,lon2)
% function [d,hdg]=t_gcdist(lat1,lon1,lat2,lon2)
% Function to calculate distance in kilometers and heading between two
% positions in latitude and longitude.
% Assumes -90 > lat > 90  and  -180 > long > 180
%    north and east are positive
% Uses law of cosines in spherical coordinates to calculate distance
% calculate conversion constants
%
%  Code from Richard Dewey.

raddeg=180/pi;
degrad=1/raddeg;
% convert latitude and longitude to radians
lat1=lat1.*degrad;
lat2=lat2.*degrad;
in1=find(lon1>180);lon1(in1)=lon1(in1)-360;
in2=find(lon2>180);lon2(in2)=lon2(in2)-360;
lon1=-lon1.*degrad;
lon2=-lon2.*degrad;
% calculate some basic functions
coslat1=cos(lat1);
sinlat1=sin(lat1);
coslat2=cos(lat2);
sinlat2=sin(lat2);
%calculate distance on unit sphere
dtmp=cos(lon1-lon2);
dtmp=sinlat1.*sinlat2 + coslat1.*coslat2.*dtmp;

% check for invalid values due to roundoff errors
in1=find(dtmp>1.0);dtmp(in1)=1.0;
in2=find(dtmp<-1.0);dtmp(in2)=-1.0;

% convert to meters for earth distance
ad = acos(dtmp);
d=(111.112) .* raddeg .* ad;

% now find heading
hdgcos = (sinlat2-sinlat1.*cos(ad))./(sin(ad).*coslat1);

% check value to be legal range
in1=find(hdgcos>1.0);hdgcos(in1)=1.0;
in2=find(hdgcos<-1.0);hdgcos(in2)=-1.0;
hdg = acos(hdgcos).*raddeg;

% if longitude is decreasing then heading is between 180 and 360
test = sin(lon2-lon1);
in1=find(test>0.0);
hdg(in1)=360-hdg(in1);
