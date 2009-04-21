function[timestr]=timefix(time);
[m n]=size(time);
timestr=datestr(time(1:m,n),'yyyymmdd.HHMM');
%t=datestr(now,'yyyymmdd');
% timestart=str2num(t);
% tend=datestr(now+3,'yyyymmdd');
% timeend=str2num(tend);%2 days for second swan run, runs from hot file

