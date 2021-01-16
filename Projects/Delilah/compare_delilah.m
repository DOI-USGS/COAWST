 compare_delilah
%
% compare delilah InWave to Observations.
%
% jcw 04Jan2021
%

% locations of the observations.
% Figure E1, E2, anbd Table E1
sites=['CM10 CM20 CM30 CM40 CM50 CM60 CM70 CM80 CM90'];
xloc=[125.06 144.99 169.97 188.94 207.41 229.25 245.00 295.21 370.12];
yloc=[985.95 985.94 985.61 985.95 985.88 986.08 985.91 985.97 986.11];
zloc=[-0.28 -1.39 -0.87 -1.32 -1.88 -2.33 -3.03 -3.68 -4.25];

%plot observation data
cd 'E:\data\models\InWave\readswan\Projects\Delilah\data\gages'

%1) DEL_Oct1990.depth	Gauge depths and local water depth
fid=fopen('DEL_Oct1990.depth');
line=fgetl(fid);
%  here is the header
%  1    2     3    4    5            6      7        8       9       10      11      12     13      14       15       16      17      18     19       20      21      22    
% Year Month Day Time   DecDay     TD231    ZC20     D20    ZC30     D30    ZC40     D40    ZC50     D50    ZC60     D60    ZC70     D70    ZC80     D80    ZC90     D90    ZC31     D31     ZC32      D32    ZC33      D33    ZC34     D34     ZC35     D35     ZC54      D54    ZC71      D71    ZC72      D72    ZC73      D73    ZC74      D74
for mm=1:100000
  line=fgetl(fid);
  if (line==-1); break; end
  depth(mm,:)=str2num(line);
end
fclose(fid);
depth(depth==-99.99)=nan;
time_depth=depth(:,5)-depth(1,5);
%  plot the data water levels
figure
hold on
plot(time_depth,depth(:, 8)+zloc(2))  % D20
plot(time_depth,depth(:,10)+zloc(3))  % D30
plot(time_depth,depth(:,12)+zloc(4))  % D40
plot(time_depth,depth(:,14)+zloc(5))  % D50
plot(time_depth,depth(:,16)+zloc(6))  % D60
plot(time_depth,depth(:,18)+zloc(7))  % D70
plot(time_depth,depth(:,20)+zloc(8))  % D80
plot(time_depth,depth(:,22)+zloc(9))  % D90
plot(time_depth,depth(:, 6),'k','linewidth',2)  % TD231 tide gage


%2) % DEL_Oct1990.wave   Wave height, period, direction (cross-shore arrray)
fid=fopen('DEL_Oct1990.wave');
line=fgetl(fid);
%  here is the header
%  1    2     3    4    5            6      7        8       9       10      11      12     13      14       15       16      17     18      19       20      21      22      23     24      25        26     27      28      29       30     31      32    33        34      35     36      37
% Year Month Day Time   DecDay       H20     P20    DR20   SNP20     H30     P30    DR30   SNP30     H40     P40    DR40   SNP40     H50     P50    DR50   SNP50     H60     P60    DR60   SNP60     H70     P70    DR70   SNP70     H80     P80    DR80   SNP80     H90     P90    DR90   SNP90
for mm=1:100000
  line=fgetl(fid);
  if (line==-1); break; end
  wave(mm,:)=str2num(line);
end
fclose(fid);
wave(wave==-99.99)=nan;
time_wave=wave(:,5)-wave(1,5);
%  plot wave heights, per, dirs
figure
for mm=1:3
  subplot(3,1,mm)
  hold on
  plot(time_wave,wave(:, 6+mm-1))  % H20 P20 DR20
  plot(time_wave,wave(:,10+mm-1))  % H30
  plot(time_wave,wave(:,14+mm-1))  % H40
  plot(time_wave,wave(:,18+mm-1))  % H50
  plot(time_wave,wave(:,22+mm-1))  % H60
  plot(time_wave,wave(:,26+mm-1))  % H70
  plot(time_wave,wave(:,30+mm-1))  % H80
  plot(time_wave,wave(:,34+mm-1))  % H90
  if mm==1; title('Wave heigths, m'); end
  if mm==2; title('Wave periods, s'); end
  if mm==3; title('Wave directions, deg'); end
end

%3) DEL_Oct1990.curUVC Cross- and longshore currents from cross-shore array
fid=fopen('DEL_Oct1990.curUVC');
line=fgetl(fid);
%  here is the header
%  1    2     3    4    5            6      7        8       9       10      11      12     13      14       15     16      17     18      19       20      21      22      23     24      25        26     27      28      29       30     31      32    33        34      35     36       37      38       39      40      41      42      43      44      45
% Year Month Day Time   DecDay      U20     V20    ZT20    SNU20   SNV20    U30     V30    ZT30    SNU30   SNV30    U40     V40    ZT40    SNU40   SNV40    U50     V50    ZT50    SNU50   SNV50    U60     V60    ZT60    SNU60   SNV60    U70     V70    ZT70    SNU70   SNV70    U80     V80    ZT80    SNU80   SNV80    U90     V90    ZT90    SNU90   SNV90
for mm=1:100000
  line=fgetl(fid);
  if (line==-1); break; end
  uvc(mm,:)=str2num(line);
end
fclose(fid);
uvc(uvc==-99.99)=nan;
time_uvc=uvc(:,5)-uvc(1,5);
figure
for mm=1:2
  subplot(2,1,mm)
  hold on
  plot(time_uvc,uvc(:, 6+mm-1))  % U20 V20
  plot(time_uvc,uvc(:,11+mm-1))  % U30
  plot(time_uvc,uvc(:,16+mm-1))  % U40
  plot(time_uvc,uvc(:,21+mm-1))  % U50
  plot(time_uvc,uvc(:,26+mm-1))  % U60
  plot(time_uvc,uvc(:,31+mm-1))  % U70
  plot(time_uvc,uvc(:,36+mm-1))  % U80
  plot(time_uvc,uvc(:,41+mm-1))  % U90
  if mm==1; title('Cross-shore current, m/s + onshore'); end
  if mm==2; title('Longshore current, m/s,pos south'); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now get model output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  locations
cd E:\data\models\InWave\readswan\delilah04
x_rho=ncread('ocean_delilah_his_00001.nc','x_rho');
y_rho=ncread('ocean_delilah_his_00001.nc','y_rho');
for mm=1:9
  xlocr(mm)=floor(interp1(x_rho(:,1),[1:1:size(x_rho,1)],xloc(mm)));
  ylocr(mm)=floor(interp1(y_rho(1,:),[1:1:size(y_rho,2)],yloc(mm)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  time
for mm=1:8
  eval(['ot',num2str(mm),'=ncread(''ocean_delilah_his_0000',num2str(mm),'.nc'',''ocean_time'');']);
end
ot=[ot1; ot2; ot3; ot4; ot5; ot6; ot7; ot8]/3600/24;
%here we offset the timess to make  equal. Their time starts at  1990  10    2  0015 275.00000
%my tiem started as 0 = sept 30 0000 so that sept 30 would be day s 00.00, then Oct 1 would be 1.xxxx
ot=ot-2;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  zeta
for jj=2:9  % number of sites
  for mm=1:8 %number of his files
    eval(['zeta',num2str(mm),'=squeeze(ncread(''ocean_delilah_his_0000',num2str(mm),'.nc'',''zeta'',[xlocr(jj) ylocr(jj) 1],[1 1 Inf]));']);
  end
  zeta(jj,:)=[zeta1; zeta2; zeta3; zeta4; zeta5; zeta6; zeta7; zeta8];
end
zeta(1,:)=nan;
plot(ot,zeta,'linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ubar
for jj=2:9  % number of sites
  for mm=1:8 %number of his files
    eval(['ubar',num2str(mm),'=squeeze(ncread(''ocean_delilah_his_0000',num2str(mm),'.nc'',''ubar'',[xlocr(jj) ylocr(jj) 1],[1 1 Inf]));']);
  end
  ubar(jj,:)=[ubar1; ubar2; ubar3; ubar4; ubar5; ubar6; ubar7; ubar8];
end
ubar(1,:)=nan;
plot(ot,ubar,'linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  vbar
for jj=2:9  % number of sites
  for mm=1:8 %number of his files
    eval(['vbar',num2str(mm),'=squeeze(ncread(''ocean_delilah_his_0000',num2str(mm),'.nc'',''vbar'',[xlocr(jj) ylocr(jj) 1],[1 1 Inf]));']);
  end
  vbar(jj,:)=[vbar1; vbar2; vbar3; vbar4; vbar5; vbar6; vbar7; vbar8];
end
vbar(1,:)=nan;
plot(ot,-vbar,'linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Hwave
for jj=2:9  % number of sites
  for mm=1:8 %number of his files
    eval(['hwave',num2str(mm),'=squeeze(ncread(''ocean_delilah_his_0000',num2str(mm),'.nc'',''Hwave'',[xlocr(jj) ylocr(jj) 1],[1 1 Inf]));']);
  end
  hwave(jj,:)=[hwave1; hwave2; hwave3; hwave4; hwave5; hwave6; hwave7; hwave8];
end
hwave(1,:)=nan;
plot(ot,hwave,'linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Dwave
for jj=2:9  % number of sites
  for mm=1:8 %number of his files
    eval(['dwave',num2str(mm),'=squeeze(ncread(''ocean_delilah_his_0000',num2str(mm),'.nc'',''Dwave'',[xlocr(jj) ylocr(jj) 1],[1 1 Inf]));']);
  end
  dwave(jj,:)=[dwave1; dwave2; dwave3; dwave4; dwave5; dwave6; dwave7; dwave8];
end
dwave(1,:)=nan;
plot(ot,90-dwave,'linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Pwave
for jj=2:9  % number of sites
  for mm=1:8 %number of his files
    eval(['pwave',num2str(mm),'=squeeze(ncread(''ocean_delilah_his_0000',num2str(mm),'.nc'',''Pwave_top'',[xlocr(jj) ylocr(jj) 1],[1 1 Inf]));']);
  end
  pwave(jj,:)=[pwave1; pwave2; pwave3; pwave4; pwave5; pwave6; pwave7; pwave8];
end
pwave(1,:)=nan;
plot(ot,90-pwave,'linewidth',2)
