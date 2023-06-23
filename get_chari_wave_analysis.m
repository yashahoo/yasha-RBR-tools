
function [wavedata hs1  hm1 S Pn Pstat]=get_chari_wave_analysis(time,data,srate,burst_duration);

% [wavedata hs1  hm1 S Pstat]=get_chari_wave_analysis(time,data,srate,burst_duration);
% burst_duration = number of seconds in procesing window (2048,4096,8192)
% time=ts;
% data=p;
% srate=2;
time1=time;
data1=data;
% PROGRAM TO INCLUDE DECIMAL TIME in RBR pressure sensors DATA

%load R077285_data.txt

%data1=R124029_20170410_1222(1:15274000,:);


[num1 b]=size(data1);
% time1=zeros(num1,1);
% 
%  % define times of busts
%   for i=1:num1
%  Month=data1(i,2);
%  Date=data1(i,3);
%  Year=data1(i,1);
%  Hour=data1(i,4);
%  Min=data1(i,5);
%  Sec=data1(i,6);
%  time1(i) = (datenum(Year,Month,Date,Hour,Min,Sec));
%   end
%  
% figure(1)
% plot(tim,data1)
% hold on 
% datetick('x',19)
% grid on
% xlabel('Time ')
% ylabel('Water Level (m)')


% wave parameters

% THE OUTPUTS ARE (columns 1-10)
%    time
%	c     - number of waves
%       Tz    - mean wave period
%	Hrms  - rms wave height
%	Hmax  - largest crest-to-trough wave
%	Tmax  - period of Hmax
%	Hs    - significant wave height; mean height of largest 1/3 of waves
%	Ts    - the mean period of the largest 1/3 of the waves
%       H10   - mean height of largest 1/10 of waves
%	T10   - the mean period of the largest 1/10 of the waves

%Pstat =
%	  [Tpeak Tot_Var Trend_Var Inf_Var Swell_Var Wind_Var Noise_Var
%	    Trend_Per Inf_Per Swell_Per Wind_Per Noise_Per
%	      Hrms Trend_ht Inf_ht Swell_ht Wind_ht Noise_ht]
%number of 2048 sections in file
% isectf=7458;
% isectf=floor(length(data)/(2048*srate/2)); % yh????? was just 2048
% isectf=floor(length(data)/(2048*srate)); % yh????? was just 2048
isectf=floor(length(data)/(burst_duration*srate)); % 

% ----------------------------------------------------------------------- %
%                       cutoff frequencies
% ----------------------------------------------------------------------- %
%isectf=10;
wavestata=zeros(30,isectf);
hs1=zeros(isectf,1);
hm1=zeros(isectf,1);

dof=4;			% required degrees of freedom for the spectra
tcut=0.005; %[0.0033]; %0.005;		% cutoff frequency for long term trend 200s [300s]
icut=0.033; %0.05;		% cutoff frequency for infragravity waves 30 s
scut=0.125; %0.15;	%8sec not 6.67 secs	% cutoff frequency for swell waves
wcut=0.5;		% cutoff frequency for wind waves
% ----------------------------------------------------------------------- %
   % starting point
   spoint = 1;

   % data sampling rate
%    srate = 2;

   % number of points
%    numpoints = 2048*srate/2;  % 
%    numpoints = 2048*srate;  % yh?????
   numpoints = burst_duration*srate;  % yh?????

   % number of points for fft
   npoints = numpoints;
   
   f = srate*(1:npoints/2)/npoints;

   % plot only for waves longer than 4 second period
   hhpoint = 0.5*npoints;

disp(['Starting with point ', num2str(spoint)])

for i=1:isectf
epoint = spoint + numpoints - 1;
disp(['Analyzing ', num2str(numpoints), ' points from ', num2str(spoint), ' to ', num2str(epoint)])

% spectral analysis and plotting
%
% % [P] = spectrum(data1(spoint:epoint,9), npoints);
% % P(:,1:2) = 2*srate/npoints*P(:,1:2);
wavestata(2:10,i)=wavepar(data1(spoint:epoint),srate,icut);


wavestata(1,i)=double(mean(time1(spoint:epoint,1)));
 [f,Pn(:,i),Pnstat]= ...
   	autospec(data1(spoint:epoint),srate,dof,tcut,icut,scut,wcut);
sigma=2*pi/wavestata(8,i);
h=mean(data1(spoint:epoint));
k= wavenumber(sigma,h);
pr=1/cosh(k*h);
wavestata(29,i)=pr;
wavestata(11:28,i)=Pnstat;
Pstat(:,i)=Pnstat;
spoint = spoint+numpoints;
end
wavedata=wavestata';

for i=1:isectf
    hs1(i)=wavedata(i,7)/wavedata(i,29);
    hm1(i)=wavedata(i,5)/wavedata(i,29);
    S(i,1:15)=datestr(wavedata(i,1),30);
       
end
% close
% 
% 
% figure(2)
% subplot(211)
% %plot(wavedata(:,1),wavedata(:,7))
% ylabel('Wave Height (m)')
% hold on
% %plot(wavedata(:,1),wavedata(:,5),'r-')
% plot(wavedata(:,1),hs1,'b-')
% plot(wavedata(:,1),hm1,'r-')
% datetick('x',19)
% axis tight
% 
% subplot(212)
% plot(wavedata(:,1),wavedata(:,8))
% datetick('x',19)
% axis tight
% ylabel('Wave Period (s)')
% xlabel('Time (Days February 2016)')
% 
% figure(3)
% 
% subplot(511)
% plot(time1(:,1),data1(:,9))
% hold on
% % axis([736657 736746 8.0 14.0])
% grid on
% xlabel('Time ')
% ylabel('Water Level (m)')
% datetick('x',19)
% axis tight
% 
% subplot(512)
% plot(wavestata(1,:),wavestata(7,:))
% hold on
% plot(wavestata(1,:),wavestata(5,:),'r-')
% xlabel('Time (days)')
% ylabel('Hs (m)')
% grid on
% % axis([736657 736746 0.0 4])
% datetick('x',19)
% axis tight
% 
% subplot(513)
% plot(wavestata(1,:),wavestata(3,:))
% hold on
% %plot(wavestata(1,:),wavestata(10,:))
% xlabel('Time (days)')
% ylabel('Tz (s)')
% grid on
% % axis([736657 736746 8 15])
% datetick('x',19)
% axis tight
% 
%  subplot(514)
%  plot(wavestata(1,:),wavestata(26,:))
%  hold on
%  plot(wavestata(1,:),wavestata(27,:),'r-')
%  xlabel('Time (days)')
%  ylabel('Swell/sea height (m)')
% grid on
% % axis([736657 736746 0.0 2])
% datetick('x',19)
% axis tight
% 
% subplot(515)
%  plot(wavestata(1,:),wavestata(20,:))
%  hold on
%  plot(wavestata(1,:),wavestata(21,:),'r-')
%  xlabel('Time (days)')
%  ylabel('Percentage (swell/sea)')
% grid on
% % axis([736657 736746 0.0 100])
% datetick('x',19)
% axis tight
% 
% 
% figure;
%  plot(wavestata(1,:),wavestata(25,:))
%  hold on
%  xlabel('Time (days)')
%  ylabel('IG height (m)')
% grid on
% % axis([736657 736746 0.0 2])
% datetick('x',19)
% axis tight
