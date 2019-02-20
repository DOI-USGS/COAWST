function [const,sat,shallow]=t_getconsts(ctime);
% T_GETCONSTS Gets constituent data structures
% [CONST,SAT,SHALLOW]=T_GETCONSTS returns data structures holding
% information for tidal analyses.
%
% Variables are loaded from 't_constituents.mat', otherwise the 
% ascii files 'tide3.dat' (provided with the IOS analysis package)
% and 't_equilib.dat' are read, and the results stored in 
% 't_constituents.mat' for future use.
%
% [...]=T_GETCONSTS(TIME) recomputes the frequencies from the 
% rates-of-change of astronomical parameters at the matlab TIME given.

% R. Pawlowicz 11/8/99
% Version 1.0


if exist('t_constituents.mat','file');
  load t_constituents
else

  nc=146;
  empvec=zeros(nc,1)+NaN;
  const=struct('name',setstr(zeros(nc,4)),... % constituent names
               'freq',empvec,...              % and frequencies (cph)
               'kmpr',setstr(zeros(nc,4)),... % names of comparisons
	       'ikmpr',empvec,'df',empvec,... % ..and their index (into .name)
	       'doodson',zeros(nc,6)+NaN,...  % doodson#s (when available)
	       'semi',empvec,...              % phase offsets
	       'isat',empvec,...              % index into "sat"
	       'nsat',empvec,...              % # of associated satellites in "sat"
	       'ishallow',empvec,...          % index in "shallow"
	       'nshallow',empvec,...          % # of generating freqs in "shallow"
	       'doodsonamp',empvec,...        % Equilibrium Amplitude (when available)
	       'doodsonspecies',empvec);      % Species
  nsat=162;
  sat=struct('deldood',zeros(nsat,3),...  % changes in last 3 doodson#s 
             'phcorr',zeros(nsat,1),...   % phase corrections
             'amprat',zeros(nsat,1),...   % amplitude corrections
	     'ilatfac',zeros(nsat,1),...  % latitude-dependent correction type
	     'iconst',zeros(nsat,1));     % index of major (in const.)
  
  nshl=251;
  shallow=struct('iconst',zeros(nshl,1),... % index of shallow constituent name
                 'coef',zeros(nshl,1),...   % corresponding combination number and
                 'iname',zeros(nshl,1));    % index of main constituent
		 
  
  
  fid=fopen('tide3.dat');
  if fid==-1,
    error('Can''t find constituent input file ''tide3.dat''!');
  end;
  l=fgetl(fid);
  k=0;
  while length(l)>24,
    k=k+1;
    const.name(k,:)=l(5:8);
    const.freq(k)=sscanf(l(14:25),'%f');
    nm=[l(30:end) '    '];
    const.kmpr(k,:)=nm(1:4); 
    l=fgetl(fid);
  end

  % Coefficients without comparison constituent are not used
  % in the present configuration.

  const.df=zeros(length(const.freq),1);
  
  for k=find(any(const.kmpr'~=' '));   
    j1=strmatch(const.kmpr(k,:),const.name);
    const.ikmpr(k)=j1;
    const.df(k)=abs(const.freq(j1)-const.freq(k));
  end
  const.df(1)=0;  % Leave df(1)=0 to remove z0 from this list

  % Skip blank lines.
  l=fgetl(fid);l=fgetl(fid);l=fgetl(fid);
 
  % Now decode the doodson# and satellite information.
  
  k=0;
  while length(l)>10,
    kon=l(7:10);
    j1=strmatch(kon,const.name);
    vals=sscanf(l(11:end),'%f');
    const.doodson(j1,:)=vals(1:6);
    const.semi(j1)=vals(7);
    if vals(8)~=0,  % Satellite data follows
      const.nsat(j1)=vals(8);
      m=vals(8);sats=[];
      while m>0,
        l=fgetl(fid);
	l=[l ' 0'];
        for n=1:min(m,3);
          if l(n*23+10)==' ',l(n*23+10)='0'; end;
	  sats=[sats,l(n*23+[-11:8]) ' ' l(n*23+10),' '];
	  m=m-1;
	end;
     end;
      vals=sscanf(sats,'%f',[6 Inf])';
      nst=size(vals,1); 
      if nst~=const.nsat(j1), error('# of satellites does not match input'); end;
      sat.deldood(k+[1:nst],:)=vals(:,1:3);
      sat.phcorr(k+[1:nst])=vals(:,4);
      sat.amprat(k+[1:nst])=vals(:,5);
      sat.ilatfac(k+[1:nst])=vals(:,6);
      sat.iconst(k+[1:nst])=j1;
      const.isat(j1)=k+1;
      k=k+nst;
    end;
    l=fgetl(fid);
  end;
  
  % Shallow water constituents - we need to get these in terms
  % of their original!
  
  l=fgetl(fid);
  k=0;
  while length(l)>3,
    kon=l(7:10);
    j1=strmatch(kon,const.name);
    nsh=sscanf(l(11:12),'%d');
    const.nshallow(j1)=nsh;
    shallow.iconst(k+[1:nsh])=j1;
    for m=1:nsh,
     shallow.coef(k+m)=sscanf(l(m*15+[0:3]),'%f');
     shallow.iname(k+m)=strmatch(l(m*15+[5:6]),const.name);
    end;
    const.ishallow(j1)=k+1;
    k=k+nsh;
    l=fgetl(fid);
  end;  
   
  %% Get the equilibrium amplitudes from Doodson's development.
  
  fid=fopen('t_equilib.dat');
  if fid==-1,
    error('Can''t find equilibrium amplitude dataset');
  end;

  % Now parse file, which is in format Name species A B.
  
  l=fgetl(fid);
  while l(1)=='%',
    l=fgetl(fid);
  end;
  while length(l)>1,
    j1=strmatch(l(1:4),const.name);
    vals=sscanf(l(5:end),'%f');
    if vals(2)~=0,
      const.doodsonamp(j1)=vals(2)/1e5;
      const.doodsonspecies(j1)=vals(1);
    else
      const.doodsonamp(j1)=vals(3)/1e5;
      const.doodsonspecies(j1)=-vals(1);
    end;
    l=fgetl(fid);
  end;
   
  save t_constituents const sat shallow
end;

if nargin==1 & ~isempty(ctime), % If no time, just take the "standard" frequencies,
                                % otherwise compute them from derivatives of astro
 [astro,ader]=t_astron(ctime);  % parameters. This is probably a real overkill - the
 ii=isfinite(const.ishallow);     % diffs are in the 10th decimal place (9th sig fig).
 const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
 for k=find(ii)',
   ik=const.ishallow(k)+[0:const.nshallow(k)-1];
   const.freq(k)=sum( const.freq(shallow.iname(ik)).*shallow.coef(ik) );
 end;
end; 



