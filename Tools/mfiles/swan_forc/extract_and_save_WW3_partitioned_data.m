%This code will read a cleaned up version of partitioned data (say for a
%range of longitudes and latitudes of your interest) and then create a SWAN
%format wave forcing file. 

%% STEP 1: Housekeeping
%clear all %#ok<*CLALL>

%% STEP 2: Get filename
%fname = ('partition_ww3_wc_200606_scb1');
%sname = 'WW3_Partitions.mat';  %This is file where we save the variables

%% STEP 3: Read data
%A     = importdata(fname);   %jcw created "A" from the extracted partition file.

%% STEP 4: Now identify all locations which have the time index (others are partitions)
if (0)
  Prt   = A(:,5);
  %j     = find((Prt-floor(Prt))==0 & Prt<10);   % This provides indices which list the number of partitions.
  j     = find(Prt==99999);
  YY_t  = A(j,1);                                % Now with the right choice of j indices we have the time.  
  ij    = YY_t<20;                               % Safety check for random values
  j(ij) = [];
  clear ij YY_t
end

%% STEP 5: Configure, time, Latitude, Longtiude, and number of partitions.
YYYYMODD = A(j,1);
HHMMSS   = A(j,2);
PtNo     = A(j,5);
Lon      = A(j,4)-360;
Lat      = A(j,3);

%% STEP 6: Identify the segments which correspond to each location (depends on number of partitions).
ID       = find(diff(HHMMSS)~=0);
SegB     = [1;[ID]+1];
SegE     = [[ID];length(HHMMSS)];
fm       = [0.02:0.01:0.5];
dth      = 10;

%% STEP 7: Loop through the code and read the partitions for each location corresponding to each time.
count = 0;

for p    = 1:1:length(SegB)
    qmax = length(Lon(SegB(p):1:SegE(p)));
    for q       = 1:1:qmax
        count   = count+1;
%       Nprt    = A(j(count),5);
        Nprt    = A(j(count),6);
%        if Nprt>0
%            Hsig    = A(j(count)+2:j(count)+1+Nprt,2);
%            Tp      = A(j(count)+2:j(count)+1+Nprt,3);
%            theta   = A(j(count)+2:j(count)+1+Nprt,5);
%            sigma   = A(j(count)+2:j(count)+1+Nprt,6);
%            gamma   = 2;
%            Hsig1   = A(j(count)+1:j(count)+1,2);
%            Tp1     = A(j(count)+1:j(count)+1,3);
%            theta1  = A(j(count)+1:j(count)+1,5);
%            sigma1  = A(j(count)+1:j(count)+1,6);
%        else
            Hsig    = A(j(count)+1:j(count)+1+Nprt,2);
            Tp      = A(j(count)+1:j(count)+1+Nprt,3);
            theta   = A(j(count)+1:j(count)+1+Nprt,5);
            sigma   = A(j(count)+1:j(count)+1+Nprt,6);
            gamma   = 2;
            Hsig1   = A(j(count)+1:j(count)+1,2);
            Tp1     = A(j(count)+1:j(count)+1,3);
            theta1  = A(j(count)+1:j(count)+1,5);
            sigma1  = A(j(count)+1:j(count)+1,6);
%        end
        [Sf_Temp,Sfth_Temp,Freq,Dir] = WWIIIpartitioned2full2D(Hsig,Tp,theta,sigma,gamma,fm,dth);
        Sfth(p,q,:,:) = Sfth_Temp;
%       str= sprintf('Doing %d and %d',p,q);
%       disp(str);
        clear Nprt Hsig Tp theta sigma gamma Sf_Temp Sfth_Temp Sfth_Temp1 Sf_Temp1
    end
end    

%% STEP 7: Fix time

YYYY = floor(YYYYMODD/10000);
MODD = YYYYMODD-YYYY*10000;
MO   = floor(MODD/100);
DD   = MODD-100*MO;
HH   = floor(HHMMSS/10000);
MMSS = HHMMSS-10000*HH;
MM   = floor(MMSS/100);
SS   = MMSS-100*MM;
Time = datenum(YYYY,MO,DD,HH,MM,SS);
Lon  = Lon(SegB(10):SegE(10));
Lat  = Lat(SegB(10):SegE(10));

%% STEP 8: Save file
eval(['save -v7.3 ',sname,' Freq Dir Sfth Time Lon Lat']);
