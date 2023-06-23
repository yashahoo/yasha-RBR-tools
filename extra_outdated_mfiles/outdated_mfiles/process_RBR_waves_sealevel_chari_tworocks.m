% yasha hetzel may 2018
% this m-file reads raw .rsk file downloaded from RBR pressure sensors
% using rsk-tools (downloadad from rbr website)
% and saves to netcdf file. then it reads the netcdf file and calculates
% waves and makes a number of plots that are saved to a directory along
% with the processed data as matfile (all) or the hourly wave data as a csv
% text file
% **There are a lot of m-files required so you need to set the path to make 
% sure all the m-files in the yasha_RBR_tools are available.
% 
% ** this was writen on a mac which uses  '/'  in the file naming, so these
% need to be changed to '\' to work on a pc

clear all;
clc;
close all;

% required:
% t-tide
% rotateXlabels
% subaxis
% exportfig
% lanczosfilter
% and some otehr including rbr rsk-tools

restoredefaultpath
% % addpath('/Users/Yasha/Documents/MATLAB/m-files')
% % addpath ( '/Users/Yasha/Documents/MATLAB/m-files/mexcdf/snctools');
% % addpath ( '/Users/Yasha/Documents/MATLAB/m-files/IVICA_matlab_tools');


% addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/subaxis'));
% addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rotateXLabels')); 
% addpath ( '/Users/Yasha/Documents/MATLAB/m-files/t_tide_v1')
% % addpath ( '/Users/Yasha/Documents/MATLAB/m-files/altmany-export_fig-76c775b')
% 
% 
% addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/yasha_RBR_tools'));
% addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rbr-rsktools'));
% addpath('t_tide_v1');
% addpath ('altmany-export_fig-76c775b')
% addpath ('subaxis')

restoredefaultpath % this just fixes if there are any confliting functions

% set the paths here relevant to your computer
addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/yashatools/yasha_RBR_tools')); % [edit]  this should contain all the functions you need, except the rsk-tooks that are required to read the raw data (see below)
% addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rbr-rsktools-2.3.0')); %  [edit] this should be installed on your computer somewhere (it has been tested ususing old version)
addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rbr-rsktools')); % use this old version to ensure it reads temperature!
%%  ------------------------settings-------------------------------------- %

% this for reading raw file
% % dir_ref='/Users/Yasha/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/SN124029_Gracetown_offshore20180519_1928/' % this is directory where raw .rsk file resides
% rawfile='/Users/Yasha/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/SN124029_Gracetown_offshore20180519_1928/SN124029_Gracetown_offshore20180519_1928.rsk' % complete path to rsk file
% % rawfile=[ 'SN124029_Gracetown_offshore20180519_1928.rsk']
location='TwoRocks'; % used for filenaming

% rawfile='\Users\Yasha\Documents\RESEARCH\DATA\MEASURED\RBRpressure_sensors\SN124029_Gracetown_offshore20180519_1928\SN124029_Gracetown_offshore20180519_1928.rsk'
% rawfile='/Users/Yasha/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/SN124029_Gracetown_boatramp_20180721/124029_20180721_1215.rsk'
% rawfile='/Users/Yasha/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/SN124028_Gracetown_Offshore_20180819_2028/124028_20180819_2028.rsk'
% rawfile='/Users/Yasha/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/SN124028_Gracetown_Offshore_20180819_2028/124028_20180819_2028.rsk'
% rawfile='/Users/Yasha/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/SN124029_Gracetown_boatramp_20181113_1053/124029_20181113_1053.rsk'
rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/SN77826_chari_TwoRocks_2019winter/077826_20200127_1413.rsk'

% deployment info
sample_rate_hz=2; % sampling in hz (normally 2)
max_period=25; % seconds for cutoff in calculating waves
% begin=datenum(2018,1,13,18,0,0); % start of good data (can get from running initial reading of data in next section)
% finish=datenum(2018,3,3,16,30,0); % end of good data (can get from running initial reading of data in next section)
begin=datenum(2019,4,28,13,0,0);
finish=datenum(2019,8,22,2,40,0);
% 28-Apr-2019 13:00:00
% 22-Aug-2019 02:40:08

% flag bad data with nans [if no bad data make badstart/end = NaN]. this
% seems to cause some problemss,so avoid for now if possible
% badstart=datenum(2018,2,4,4,17,3);
%   badend=datenum(2018,2,4,4,19,20);
badstart=nan;
badend=nan;

% set where to save figures and data
fs=12; % set font size for figs
whereput=[rawfile(1:end-4) '_processed'];
mkdir(whereput);

% options
% plot weekly highpass filter (was useful to look at boatwakes in albany)
plot_highpass_weeks='N'  % 'N' if don't want to plot

% save the processed data to matfile
save_matfile='N'

% read raw file and save to netcdf

% 
% 
[ncref]=read_RBRyh(rawfile,location,whereput);

% % % this used later
% % ncref=[whereput '/' rawfile(1:end-4) '.nc'];
% % ncref=[whereput '/' rawfile(end-23:end-4) '.nc'];
% ncref=[whereput '/' rawfile(end-23:end-4) '.nc'];
close all;
% /Volumes/Margs_Clone/RBRpressure_sensors/SN77826_chari_TwoRocks_2019winter/077826_20200127_1413_processed/077826_20200127_1413.nc
% /Volumes/Margs_Clone/RBRpressure_sensors/SN77826_chari_TwoRocks_2019winter/077826_20200127_1413_processed/077826_20200127_1413.nc


%% read the nc file and begin processing
% ts=ncread([dir_ref ncref],'time');%+datenum(2000,1,1);
% p=ncread([dir_ref ncref],'press');
ts=ncread(ncref,'time');%+datenum(2000,1,1);
p=ncread(ncref,'press');
instrument=ncreadatt(ncref,'/','Serial_Number'); % serial number here for id

try
% temp=ncread([dir_ref ncref],'temperature');
temp=ncread([ncref],'temperature');
catch
    disp('Temperature not available')
end

% trim to good data
% figure; 
% plot(ts,p,'r')

% interactively select start and finish times
if isnan(begin) && isnan(finish)
[begin finish]=get_beginfinish(ts, p);
end


badi=find(ts<begin | ts>finish); %SN124096

% interpolate across small data gaps
if ~isnan(badstart);
badp=find(ts>badstart & ts<badend);
p(badp)=nan;
p=interp1gap(p);
end

ts(badi)=[];
p(badi)=[];

if exist('temp','var')
temp(badi)=[];
end

figure;
hold on
plot(ts,p,'b')
datetick
datetick('x','yyyy-mm-dd HH:MM:SS','keeplimits')

% close
%% calculate waves
clear Hs Ts Htime Hrms Hmax Tmax H10 T10 Tz


sph=sample_rate_hz*3600;




numdays=floor(ts(end)-ts(1));
clear waves
co=0;
for hr=1:numdays*24
   co=co+1; 
   inds=[hr*sph-(sph-1):hr*sph];
psub=p(inds);
tsub=ts(inds);


clear waves
waves=wavepar_yh(psub,sample_rate_hz,1/max_period); %2 hz

Hs(co)=waves.Hs;
Ts(co)=waves.Ts;
Htime(co)=tsub(sph/2+1);
% datestr(tsub);

Hpress(co)=nanmean(psub); % hourly pressure

Hrms(co)=waves.Hrms;
Hmax(co)=waves.Hmax;
Tmax(co)=waves.Tmax;
H10(co)=waves.H10;
T10(co)=waves.T10;
Tz(co)=waves.Tz;

if exist('temp','var')
    tempsub=nanmean(temp(inds));
Htemperature(co)=tempsub;
end

end

xlims=[ts(1) ts(end)]; % all
xlims=[begin finish];% good
xt=[xlims(1):xlims(2)];
yt=[0:2:30];





f=figure; 
f.Position=[66 207 1600 450];
f.Color='w';
% set(gcf,'Position',[200 100 700 300])
% set(gcf,'PaperPositionMode','auto')
plot(Htime,Hs)
hold on
plot(Htime,H10)
plot(Htime,Hmax)
legend('Hs','H10','Hmax')
ylabel('Wave height (m)')
set(gca,'xlim',xlims,'xtick',xt)
% set(gca,'ylim',[0 ceil(max(Hmax))])
set(gca,'ylim',[0 .1+max(Hmax)])
datetick('x','mm/dd','keeplimits','keepticks')
rotateXLabels( gca(), 90 )
title(instrument)
fname=[whereput '/' instrument '_WAVE_HEIGHT_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
% % print('-dpng','-r200', fname);
% print('-dpdf', '-painters', fname);
% % print('-depsc', '-painters', fname);
export_fig(fname,'-pdf','-png','-transparent')

% figure; 
% set(gcf,'Position',[200 100 700 300])
% set(gcf,'PaperPositionMode','auto')
f=figure; 
f.Position=[66 207 1600 450];
f.Color='w';
plot(Htime,T10)
hold on
plot(Htime,Ts)
set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
plot(Htime,Tz)
legend('T10','Ts','Tz')
ylabel('Period (s)')
datetick('x','mm/dd','keeplimits','keepticks')
rotateXLabels( gca(), 90 )
fname=[whereput '/' instrument '_WAVE_PERIOD_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
% print('-dpng','-r200', fname);
% print('-dpdf', '-painters', fname);

 export_fig(fname,'-pdf','-png','-transparent')
 
 
 if exist('temp','var')
%      figure;
%      set(gcf,'Position',[200 100 700 300])
%      set(gcf,'PaperPositionMode','auto')
    f=figure;
    f.Position=[66 207 1600 450];
    f.Color='w';
    plot(Htime,Htemperature)
    hold on
%      legend('Hs','H10','Hmax')
     ylabel('temperature ( deg C)')
     set(gca,'xlim',xlims,'xtick',xt)
     set(gca,'ylim',[floor(min(Htemperature)) ceil(max(Htemperature))])
     datetick('x','mm/dd','keeplimits','keepticks')
     rotateXLabels( gca(), 90 )
     title(instrument)
     fname=[whereput '/' instrument '_TEMPERATURE_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
     % % print('-dpng','-r200', fname);
     % print('-dpdf', '-painters', fname);
     % % print('-depsc', '-painters', fname);
     export_fig(fname,'-pdf','-png','-transparent')
 end
%% get IG waves (this overlaps with above)

[wavedata hs1  hm1 S]=get_chari_wave_analysis(ts,p,2);

% figure; 
% set(gcf,'Position',[200 100 700 300])
% set(gcf,'PaperPositionMode','auto')
f=figure;
f.Position=[66 207 1600 450];
f.Color='w';
hold on
box on
plot(wavedata(:,1),wavedata(:,25))
set(gca,'xlim',xlims,'xtick',xt)
 ylabel('IG height (m)')
datetick('x','mm/dd','keeplimits','keepticks')
rotateXLabels( gca(), 90 )

% figure;
%  plot(wavestata(1,:),wavestata(25,:))
%  hold on
%  xlabel('Time (days)')
%  ylabel('IG height (m)')
fname=[whereput '/' instrument '_IG_WAVES_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
% print('-dpng','-r200', fname);
% print('-dpdf', '-painters', fname);

 export_fig(fname,'-pdf','-png','-transparent')
 
 
 
 % get hourly ig
 tig=wavedata(:,1);
 ig=wavedata(:,25);
 
 Hig=interp1(tig, ig,Htime);
 
 
%% get lowpass sea level (lanczos or butterworth 38 hours)

% SN124028
% p(1:289286)=0;
% p(9325510:end)=0;

% pad with zeros to make sure start is ok
Hpress0=Hpress;
% Hpress0(1:8)=0;
% Hpress0(1294:end)=0;
Hpress0=[zeros(1,100) Hpress0 zeros(1,100)]; 

low_cutoff=38;
% lowpass filter model
[lp_wl,B,A]=lp_filter(6,low_cutoff,1,Hpress0);
Hpress0(1:100)=[];Hpress0(end-99:end)=[];
lp_wl(1:100)=[];lp_wl(end-99:end)=[];
lp_wl(1:28)=nan;
lp_wl(1029:end)=nan;
Hpress0(1:28)=NaN;
Hpress0(1029:end)=NaN;

figure; 
hold on;
plot(Htime,Hpress0)
plot(Htime,lp_wl)
datetick('x','yyyy/mm/dd HH:MM:SS','keeplimits')

figure; 
hold on;
plot(Htime,Hpress-lp_wl)
% plot(Htime,lp_wl)
datetick('x','yyyy/mm/dd HH:MM:SS','keeplimits')


% [lp_wl,B,A]=lp_filter(6,3600,1,p(289286:289286+24*7200));

% try lanczos fulter
      t=Htime;
      xn=Hpress;
      dT=1;
      Tc = 38; % hr
  [lanczos_wl,c,h,Cx,f] = lanczosfilter(xn,dT,1/Tc,[],'low');

%   lanczos_wl=xs;
%   
  
% choose which to use (they are the same from both filters)
lowpass=lanczos_wl-nanmean(lanczos_wl);
% badi=[];
% badi=find(Htime<datenum(2017,9,22,0,0,0) | Htime>datenum(2017,11,11,11,0,0)); %SN124096
% lowpass(badi)=nan;
% lowpass=lowpass-nanmean(lowpass);
%% get bandpass IG badns ????
% 
% t_step=1/2
% cut1=0.002 
% cut2=0.45
% red=2
% % x=p-nanmean(p);
% % x(isnan(x))=0;
% 
% f_nyquist=1/t_step/2;
% f1_cutof=1/cut1/f_nyquist;
% f2_cutof=1/cut2/f_nyquist;
% % disp(['Cut off frekvencija je :'  num2str(f_cutof)]);
% % [B,A]=butter(red,[f1_cutof f2_cutof],'bandpass');
% [B,A]=butter(red,[f1_cutof f2_cutof],'stop');
% 
% 
% freqz(B,A)
% 
% y=filtfilt(B,A,x);
% 
% %
% 
% %% junk
% 
% Fsp = 2;                                     % Create Data
% Fn = Fsp/2;
% Wp = [1.0 49]/Fn;
% Ws = [0.5 50]/Fn;
% Rp=10;
% Rs=30; 
% [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs); 
% [z,p,k] = cheby2(n,Rs,Ws);
% [sos,g] = zp2sos(z,p,k);
% figure(1)
% freqz(sos, 2^16, Fsp)
% set(subplot(2,1,1), 'XLim',[0 100])             % ?Zoom? X-Axis
% set(subplot(2,1,2), 'XLim',[0 100])             % ?Zoom? X-Axis


%%

% %
% % [y,B,A]=bandpass_filter(red,cut1,cut2,t_step,x)
% 
% samplehertz=sample_rate_hz;
% sampleperiod=1/samplehertz;
% 
% % psub=p;
% % tsub=ts;
% % wt=ts;
% 
% data=p;
% data=data-mean(data);		% remove the mean
% % data=highpass(data,2,1/max_period);	% highpass filter
% 
%  [y,B,A]=bandpass_filter(5,30,300,sampleperiod,data);
% 
% [doy,frac] = date2doy(datenum(ts));
% td=datetime(datevec(wt),'format','eee yyyy-MM-dd HH:mm');
% tf = isweekend(td);
% wz=zeros(length(wt),1);
% 
% % plot weeks at a time
% for i=[doy(1):7:doy(1)+floor(doy(end)-doy(1))];
%     temp=i+datenum(2017,1,1);
% inds=find(ts>temp & ts<temp+6);
% i
% figure; 
% % set(gcf,'Position',[200 100 1200 800])
% set(gcf,'Position',[200 100 1200 400])
% set(gcf,'PaperPositionMode','auto')
% hold on
% plot(tsub,data);
% plot(wt(tf),wz(tf)-.3,'r.')
% set(gca,'xlim',[ts(inds(1)) ts(inds(end))]);
% set(gca,'xticklabel',[]);
% datetick('x','mm/dd HH:MM:SS','keeplimits')
% % rotateXLabels( gca(), 90 )
% 
% ylabel('Height (m)')
% title([instrument ' highpass (<' num2str(max_period) ' s) sea level'])
% 
% set(gca,'fontsize',12)
% 
% 
% % pause
% 
% fname=[whereput '/' instrument '_HIGHPASS_weekly_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
% % export_fig(fname,'-pdf','-png','-transparent')
% export_fig(fname,'-png','-transparent')
% close
% 
% end
% 
%  
% 
% figure; 
% hold on
% % set(gcf,'Position',[200 100 1200 800])
% set(gcf,'Position',[200 100 1200 400])
% set(gcf,'PaperPositionMode','auto')
% plot(tsub,data);
% % set(gca,'xlim',[ts(1) ts(end)]);
% % set(gca,'xticklabel',[]);
% % datetick('x','yyyy/mm/dd ','keeplimits')
% plot(wt(tf),wz(tf)-.4,'r.')
% set(gca,'xlim',xlims,'xtick',xt)
% set(gca,'ylim',[-.4 .4],'ytick',[-.6:.1:.6])
% datetick('x','mm/dd','keeplimits','keepticks')
% rotateXLabels( gca(), 90 )
% ylabel('Height (m)')
% title([instrument ' highpass (<' num2str(max_period) ' s) sea level'])
% set(gca,'fontsize',12)
% box on
% grid on
% 
% fname=[whereput '/' instrument '_HIGHPASS_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
% export_fig(fname,'-png','-transparent')




%% get tidal water levels by subtracting lowpass

% Htide=Hpress-lp_wl;

%% or do t-tide to get tides

     % predict  tide & get residual from model data
        clear names freq tidecon Htide

[names,freq,tidecon,Htide]=t_tide(Hpress,'interval',1,'start',Htime(1),'synthesis',1);
% mYout=t_predic(mts_sub,mnames,mfreq,mtidecon,'synthesis',1);
    
% calculate non-tidal residual
residual=Hpress-nanmean(Hpress)-Htide;


% figure; 
% set(gcf,'Position',[200 100 700 300])
% set(gcf,'PaperPositionMode','auto')
    f=figure;
    f.Position=[66 207 1600 450];
    f.Color='w';
hold on
plot(Htime,Hpress-nanmean(Hpress),'k');
plot(Htime,Htide)
plot(Htime,residual,'r')
% plot(Htime,lp_wl-nanmean(lp_wl),'g')
plot(Htime,lanczos_wl-nanmean(lanczos_wl),'g');
% legend('Total  water level','Tide')
ylabel('Sea level (m)')
legend('Total','Tide','Residual','lowpass') % ,'lanczos lowpass'
set(gca,'xlim',xlims,'xtick',xt)
datetick('x','mm/dd','keeplimits','keepticks')
rotateXLabels( gca(), 90 )
box on;
grid on;
fname=[whereput '/' instrument '_SEALEVEL_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
% print('-dpng','-r200', fname);
% print('-dpdf', '-painters', fname);

 export_fig(fname,'-pdf','-png','-transparent')
 
 

%% plot all

% set  date ranges
xlim(1,:)=[begin finish]; % all
% xlim(2,:)=[datenum(2017,9,4)  datenum(2017,9,10)];
% xlim(3,:)=[datenum(2017,9,13)  datenum(2017,9,16)];
% xlim(4,:)=[datenum(2017,10,1)  datenum(2017,10,25)];
% xlim(5,:)=[datenum(2017,10,13)  datenum(2017,10,16)];

clear xlims

for i=1 %:5
xlims=xlim(i,:);

xt=[xlims(1):xlims(2)];

figure;
set(gcf,'Position',[200 100 600 800])
set(gcf,'PaperPositionMode','auto')
subaxis(5,1,1)
plot(Htime,Hmax)
hold on;
plot(Htime,Hs)
plot(Htime,Hrms)
legend('Hmax','Hs','Hrms')
set(gca,'xlim',xlims,'xtick',xt)
set(gca,'ylim',[0 .1+max(Hmax)],'ytick',[0:.1:5]) %[0:.5:5] 'ylim',[0 ceil(max(Hmax))]
set(gca,'xticklabel',[])
title(instrument)
ylabel('Height (m)')

subaxis(5,1,2)
plot(Htime,Ts)
hold on;
plot(Htime,T10)
plot(Htime,Tz)
legend('Tz','T10','Ts') %'Hmax','Hs','Tz') %,'location','best'
set(gca,'xlim',xlims,'xtick',xt)
set(gca,'ylim',[5 20],'ytick',[0:2:30]) %[4 22]
set(gca,'xticklabel',[])
ylabel('Period (s)')

subaxis(5,1,3)
hold on;

plot(Htime,Hpress-nanmean(Hpress))
plot(Htime,lowpass)
plot(Htime,residual,'r')
plot([Htime(1) Htime(end)],[0 0],'color',[.5 .5 .5],'linewidth',0.5)
legend('sea level','lowpass sl','residual') % ,'location','best'
set(gca,'xlim',xlims,'xtick',xt)
set(gca,'ylim',[-.5 .5],'ytick',[-.5:.25:.5])
set(gca,'xticklabel',[])
ylabel('Height (m)')

subaxis(5,1,4)
hold on;
% plot(Htime,Hpress-lp_wl)

plot(Htime,Htide)
plot([Htime(1) Htime(end)],[0 0],'color',[.5 .5 .5],'linewidth',0.5)
legend('Tide')
set(gca,'xlim',xlims,'xtick',xt)
set(gca,'ylim',[-.5 .5],'ytick',[-.5:.25:.5])
% datetick('x','mm/dd','keeplimits','keepticks')
set(gca,'xticklabel',[])
% rotateXLabels( gca(), 90 )
ylabel('Height (m)')

subaxis(5,1,5)
hold on;
plot(wavedata(:,1),wavedata(:,25))
plot([Htime(1) Htime(end)],[0 0],'color',[.5 .5 .5],'linewidth',0.5)
legend('IG waves','location','best')
set(gca,'xlim',xlims,'xtick',xt)
set(gca,'ytick',[-1:.2:1],'ylim',[0 1])
datetick('x','mm/dd','keeplimits','keepticks')
rotateXLabels( gca(), 90 )
ylabel('Height (m)')




fname=[whereput '/' instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
% print('-dpng','-r200', fname);
% print('-dpdf', '-painters', fname);
% print('-depsc', '-painters', fname);
 export_fig(fname,'-pdf','-png','-transparent')

end 

%% save to text

Hs=Hs';Sealevel=Hpress'; Hrms=Hrms';H10=H10';Hmax=Hmax'; Hlowpass_wl=lp_wl'; 
Htide=Htide'; Hresidual=residual';
Ts=Ts';Tz=Tz';T10=T10';Tmax=Tmax'; Hig=Hig';

if exist('Htemperature','var')
Htemperature=Htemperature';
end
% varnames={'Time WST','Hrms',

% table
timestr=datestr(Htime','yyyy-mm-dd HH:MM:SS');
Time_WST=cellstr(timestr);
if exist('Htemperature','var')
T = table(timestr,Sealevel,Hlowpass_wl,Htide,Hresidual,Hs,Hrms,H10,Hmax,Hig,Ts,Tz,T10,Tmax,Htemperature);
else
T = table(timestr,Sealevel,Hlowpass_wl,Htide,Hresidual,Hs,Hrms,H10,Hmax,Hig,Ts,Tz,T10,Tmax);    
end
% % outf_tname=[whereput '/' instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '.csv'];
% outf_tname=[dir_ref instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '.csv'];
outf_tname=[whereput '/' instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '.csv'];
writetable(T,outf_tname,'Delimiter',',','QuoteStrings',true)


%% get spectra


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the spectral analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data=p(goodi); % for depth data (observed)
data=p;

% data=zeros(2431001,1);
% data=zeros(length(P),1);

% data=AnDepthmm/1000.; % for conversion frommm to meters depth

% sample frequency 5 min = 300 sec
% sampleperiod=5*60;
% sampleperiod=3600;
samplehertz=sample_rate_hz;
sampleperiod=1/samplehertz;

% loop to find first power of 2 greater than length(loc)
% the length of the north and east files should be the same
i=1;
dir=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% observations
j=1;
if dir==1
    while 2^j < length(data);
    j=j+1;
    end
else
    while 2^j < length(data);
    j=j+1;
    end
end

n2=2^j;

% to produce spectral densities and frequencies
% using a cosine taper window
% with a 95% confidence interval

out=oppsd(sampleperiod,data,1,1,64,1.095,1,1);       
f=out(:,1);
Sxx = out(:,2);

 %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do plots of spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To format the plot - Sxx vs f (log-log)

% figure;
    f=figure;
    f.Position=[66 207 900 900];
    f.Color='w';
%subplot(221)
loglog(f,Sxx,'k-','linewidth',1.5)
hold on
%axis([1e-7 1e-2 1e0 1e4])
ylabel('Spectral Density (m^2/Hz)')
xlabel('Frequency (Hz)')
% set(gca,'xlim',[min(f) max(f)]);
set(gca,'xlim',[min(f) max(f)]);
% title(['Tide gauge spectra for: ' SID2{fi} ' , ' datestr(ts(1),'dd-mm-yyyy') ' to ' datestr(ts(end),'dd-mm-yyyy') ])
ylim=get(gca,'ylim');
ymiddle=1e-4;
ytop=1e3;

% to add a vertical line at a period of 24 hours (f24)

f24=1/(24*60*60);
plot([f24 f24],ylim,'k:')
text(f24,ymiddle,' 24h','fontsize',[12])

% to add a vertical line at a period of 12 hours (f12)

f12=1/(12*60*60);
plot([f12 f12],ylim,'k:')
text(f12,ymiddle,' 12h','fontsize',[12])

% to add a vertical line at a period of 3 hours (f3)
% 
% 10 mins
f3=1/((10/60)*60*60);
plot([f3 f3],ylim,'k:')
text(f3,ytop,'10 min','fontsize',[12]) 

% f3=1/((1/60)*60*60);
% plot([f3 f3],ylim,'k:')
% text(f3,ytop,'1 min','fontsize',[12])

f20=1/300;
plot([f20 f20],ylim,'k:')
text(f20,100,'300 sec','fontsize',[12])
% 
f20=1/30;
plot([f20 f20],ylim,'k:')
text(f20,1,'30 sec','fontsize',[12])
% 
f12=1/15;
plot([f12 f12],ylim,'k:')
text(f12,.1,'15 sec','fontsize',[12])

fname=[whereput '/' instrument '_SPECTRA_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
 export_fig(fname,'-pdf','-png','-transparent')

 %% do spectrogram
% % 
% % % X=data(1:floor(length(data)/10));
% % X=p(1:length(p)); %(1:1024*100); %data; %(1:1024*100)
% % % X=data;
% % 
% % % % % % S = spectrogram(X,2*3600*38,floor(2*3600*38*.20)) ;
% % % % % % j=1;
% % % % % % while 2^j < length(X);
% % % % % %  j=j+1;
% % % % % % end
% % % % % % 
% % % % % % n2=2^j;
% % % % % 
% % % % % % to produce spectral densities and frequencies
% % % % % % using a cosine taper window
% % % % % % with a 95% confidence interval
% % % % % 
% % % % % out=oppsd(sampleperiod,data,1,1,64,1.095,1,1);
% % windlength=4096 %1*3600*6 %256 %1024*40; %
% % window=hann(windlength);
% % overlap=round(windlength*.75);  % too much over lap bad
% % [s,w,t]=spectrogram(X,window,overlap,[],2,'yaxis');
% %  tt=ts(round(t));
% % 
% %  A = abs(s);
% %  phi = angle(s);
% %    
% %    figure; 
% %    pcolor(t,log10(w),log10(A));
% % %     pcolor(tt,w,log10(A));
% %    shading flat
% %    caxis([-1 3])
% %    colormap(jet)
% %    colorbar
% % %    set(gca,'xlim',xlims,'xtick',xt)
% % % datetick('x','mm/dd','keeplimits','keepticks')
% % %    rotateXLabels( gca(), 90 )
% %  
% % 
% % 
% % %    figure; 
% % %    pcolor(t,log(w),log(A));
% % % %     pcolor(tt,w,log10(A));
% % %    shading interp
% % %    caxis([-1 3])
% % %    colormap(jet)
% % %    colorbar

%% use chari's function to calculate time-frequency plot


% [logf logj]=get_tidespecft(time,data,numpoints,overlap,samplehertz);

[logf logj timef]=get_tidespecft(ts,p,4096,1024,2);

% plot it


figure;
set(gcf,'Position',[200 100 900 500])
set(gcf,'PaperPositionMode','auto')
pcolor(timef,logf,logj(:,1:length(timef)))
shading('interp')

hold on

%axis([735600 735965  -12 -6])
% to add a vertical line at a period of 24 hours (f24)
ymiddle=0.20;
xlim=get(gca,'xlim');

f24=log(1/(12*60*60));
plot(xlim,[f24 f24],'k-')
text(timef(100),f24+0.25,' 12h','fontsize',[12])

f24=log(1/(24*60*60));
plot(xlim,[f24 f24],'k-')
text(timef(100),f24+0.25,' 24h','fontsize',[12])

% f24=log(1/(6*60*60));
% %plot(xlim,[f24 f24],'k-')
% text(timef(100),f24+0.25,'6hr','fontsize',[12])

% f24=log(1/2);
% plot(xlim,[f24 f24],'k-')
% text(timef(100),f24+0.25,' 2s','fontsize',[12])

% f24=log(1/4.5);
% plot(xlim,[f24 f24],'k-')
% text(timef(100),f24+0.25,' 4.5s','fontsize',[12])

f24=log(1/8);
plot(xlim,[f24 f24],'k-')
text(timef(100),f24+0.25,' 8s','fontsize',[12])

f24=log(1/15);
plot(xlim,[f24 f24],'k-')
text(timef(100),f24+0.25,' 15s','fontsize',[12])

f24=log(1/30);
plot(xlim,[f24 f24],'k-')
text(timef(100),f24+0.25,' 30s','fontsize',[12])

% f24=log(1/130);
% plot(xlim,[f24 f24],'k-')
% text(timef(1000),f24+0.25,'130s','fontsize',[12])

% f24=log(1/200);
% plot(xlim,[f24 f24],'k-')
% text(timef(1000),f24+0.25,'200s','fontsize',[12])


f24=log(1/300);
plot(xlim,[f24 f24],'k-')
text(timef(100),f24+0.25,' 300s','fontsize',[12])
% 
% f24=log(1/600);
% plot(xlim,[f24 f24],'k-')
% text(timef(100),f24+0.25,' 600s','fontsize',[12])

xlabel('Time')
ylabel('Frequency(Hz)')

colormap(jet)
set(gca,'xlim',[timef(1) timef(end)],'xtick',xt,'fontsize',fs)
datetick('x','mm/dd','keeplimits','keepticks')
rotateXLabels( gca(), 90 )
% set(gca,'xlim',[timef(1) timef(end)],'fontsize',fs);



fname=[whereput '/' instrument '_TIME_FREQ_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
 export_fig(fname,'-png','-transparent')


% 'VerticalAlignment', 'top'

 %% %% do highpass filter



% inds=1:length(p);
% psub=p(inds);
% tsub=ts(inds);
% inds=find(ts>begin & ts<finish); %27


psub=p;
tsub=ts;
wt=ts;

% waves=wavepar_yh(psub,2,1/10);
data=psub;
data=data-mean(data);		% remove the mean
data=highpass(data,2,1/max_period);	% highpass filter

[doy,frac] = date2doy(datenum(ts));
td=datetime(datevec(wt),'format','eee yyyy-MM-dd HH:mm');
tf = isweekend(td);
wz=zeros(length(wt),1);

% plot weeks at a time

switch plot_highpass_weeks
 case 'Y'
for i=[doy(1):7:doy(1)+floor(doy(end)-doy(1))];
    temp=i+datenum(2018,1,1);
inds=find(ts>temp & ts<temp+6);
i
figure; 
% set(gcf,'Position',[200 100 1200 800])
set(gcf,'Position',[200 100 1200 400])
set(gcf,'PaperPositionMode','auto')
hold on
plot(tsub,data);
plot(wt(tf),wz(tf)-.12,'r.')
set(gca,'xlim',[ts(inds(1)) ts(inds(end))],'ylim',[floor(min(data)) ceil(max(data))]);
set(gca,'xticklabel',[]);
datetick('x','mm/dd HH:MM:SS','keeplimits')
% rotateXLabels( gca(), 90 )

ylabel('Height (m)')
title([instrument ' highpass (<' num2str(max_period) ' s) sea level'])

set(gca,'fontsize',12)


% pause

fname=[whereput '/' instrument '_HIGHPASS_weekly_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
% export_fig(fname,'-pdf','-png','-transparent')
export_fig(fname,'-png','-transparent')
close

end
end
 %

f=figure; 
hold on
% set(gcf,'Position',[200 100 1200 800])
f.Position=[200 100 1200 400];
f.Color='w';
plot(tsub,data);
% set(gca,'xlim',[ts(1) ts(end)]);
% set(gca,'xticklabel',[]);
% datetick('x','yyyy/mm/dd ','keeplimits')
plot(wt(tf),wz(tf)-1,'r.')
set(gca,'xlim',xlims,'xtick',xt)
% set(gca,'ylim',[floor(min(data)) ceil(max(data))],'ytick',[-1:.25:2])
set(gca,'ylim',[min(data) max(data)],'ytick',[-1:.25:2])
datetick('x','mm/dd','keeplimits','keepticks')
rotateXLabels( gca(), 90 )
ylabel('Height (m)')
title([instrument ' highpass (<' num2str(max_period) ' s) sea level'])
set(gca,'fontsize',12)
box on
grid on

fname=[whereput '/' instrument '_HIGHPASS_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
export_fig(fname,'-png','-transparent')


%% save all data to mat file

switch save_matfile
    case 'Y'
% fname=[dir_ref instrument '_PROCESSED_DATA_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];

fname=[whereput '/' instrument '_PROCESSED_DATA_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
save(fname)
disp(['All processed variables saved to: ' fname '.mat'])

end

%% get time-spectra
% clear F SXX
% for hr=1:2 *24 %numdays*24
%    co=co+1; 
%    inds=[co*24*3600-(24*3600)+1:co*3600*24];
% psub=p(inds);
% tsub=ts(inds);
% 
% 
% % data=zeros(2431001,1);
% % data=zeros(length(P),1);
% 
% % data=AnDepthmm/1000.; % for conversion frommm to meters depth
% 
% % sample frequency 5 min = 300 sec
% % sampleperiod=5*60;
% % sampleperiod=3600;
% samplehertz=sample_rate_hz;
% sampleperiod=1/samplehertz;
% 
% % loop to find first power of 2 greater than length(loc)
% % the length of the north and east files should be the same
% i=1;
% dir=1;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % observations
% j=1;
% if dir==1
%     while 2^j < length(data);
%     j=j+1;
%     end
% else
%     while 2^j < length(data);
%     j=j+1;
%     end
% end
% 
% n2=2^j;
% 
% % to produce spectral densities and frequencies
% % using a cosine taper window
% % with a 95% confidence interval
% 
% out=oppsd(sampleperiod,data,1,1,64,1.095,1,1);       
% f=out(:,1);
% Sxx = out(:,2);
% 
% 
% F(:,co)=f;
% SXX(:,co)=Sxx;
% 
% 
% end
% %% plot
% figure; 
% pcolor([1:size(SXX,2)],f,log10(SXX)); 
% shading flat
% colormap(jet);
% caxis([-8 -3])
