% process DOT 5 min tide gauge data
% yasha hetzel
% 20221020


clear all; close all;clc

%% settings  

infile='/Users/00068592/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/20220818_RBR_ABB_DUNS_EB_StuBarr/GEO_Jan-Aug_2022.txt';

startdate=datenum(2022,1,1);
finishdate=datenum(2022,9,22);
whereput='/Volumes/Margs_Clone/RBRpressure_sensors/202201-09_PortGeo_tidegauge/'; %202207-09_PortGeo_tidegauge
mkdir(whereput);

save_matfile       = 'N'            % save the processed data to matfile
save_matfig        = 'N';           % save the .fig files?
print_fig          = 'Y';           % print to image file?
use_export_fig     = 'N'            % if print_fig ='Y' -->  'N' uses native matlab print function (default) OR 'Y' uses export_fig to make nice plots
img_type           = {'pdf','png'}; % if print_fig ='Y' -->  {'pdf','png','eps'} % options. can be multiple. pdf will not save if native matlab unless edit function

%% read data

[z, t, junk] = importDOTfile(infile, [52, Inf]);
z(isnan(z))=[];t(isnan(t))=[];
tt=num2str(t,'%8.4f');
mtime=datenum(tt,'yyyymmdd.HHMM');

z=z./100;
%% plot data


figure; 
plot(mtime,z);
datetick;


% 20220701-20220816 rbr deployments
% keepi=find(mtime>datenum(2022,7,1) &  mtime<datenum(2022,8,17));
% keepi=find(mtime>datenum(2022,5,1) &  mtime<datenum(2022,6,1)); % may storm
keepi=find(mtime>startdate &  mtime<finishdate);

t2=mtime(keepi);
z2=z(keepi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    do t-tide to get tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(genpath('t_tide_v1'))
addpath(genpath(pwd))
% predict  tide & get residual from model data
clear names freq tidecon tide

[names,freq,tidecon,tide]=t_tide(z,'interval',5/60,'start',mtime(1),'synthesis',1);
% mYout=t_predic(mts_sub,mnames,mfreq,mtidecon,'synthesis',1);
% calculate non-tidal residual
residual=z-nanmean(z)-tide;

% % get higher frequency tide output matching subsampled
% % time_subsampled press_subsampled
% tide_subsampled=t_predic(time_subsampled,names,freq,tidecon,'synthesis',1);
% % calculate non-tidal residual
% residual_subsampled=press_subsampled-nanmean(press_subsampled)-tide_subsampled;


[names2,freq2,tidecon2,tide2]=t_tide(z2,'interval',5/60,'start',t2(1),'synthesis',1);
% mYout=t_predic(mts_sub,mnames,mfreq,mtidecon,'synthesis',1);
% calculate non-tidal residual
residual2=z2-nanmean(z2)-tide2;

%%  lowpass filter


dT=5/60;
Tc = 1; % hr
[lanczos_wl,c,h,Cx,f] = lanczosfilter(z2,dT,1/Tc,[],'low');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create Hourly averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Create hourly averages...')

vec=datevec(t2(1));
v5  = vec(1,5)+vec(1,6)/60;
vec(5) = ceil(v5/30)*30;
vec(6) = 0;
Hstart=datenum(vec)+30/60/24;% this ensures correct number of samples at start

vec=datevec(t2(end));
v5  = vec(1,5)+vec(1,6)/60;
vec(5) = floor(v5/30)*30;
vec(6) = 0;
Hend=datenum(vec)-30/60/24; % this ensures correct number of samples at  end

Htime=[datenum(Hstart):1/24:datenum(Hend)];

Hz=[]; clear Htemperature
for i=1:length(Htime);
    idx=find(t2>Htime(i)-30/60/24 & t2<Htime(i)+30/60/24);
    Hz(i)=mean(z2(idx));

end


[Hnames,Hfreq,Htidecon,Htide]=t_tide(Hz,'interval',1,'start',Htime(1),'synthesis',1);
% mYout=t_predic(mts_sub,mnames,mfreq,mtidecon,'synthesis',1);
% calculate non-tidal residual
Hresidual=Hz-nanmean(Hz)-Htide;


%% Plots


% figure; plot(mtime,tide);datetick
% figure; plot(mtime,residual);datetick

f=figure;
f.Color='w';
set(gcf,'Position',[200 100 1600 640])
set(gcf,'PaperPositionMode','auto')
hold on

plot(t2,z2)
plot(t2,tide2+mean(z2));
plot(t2,lanczos_wl);
plot(Htime,Hz);
plot(t2,residual2,'r');
plot(Htime,Hresidual,'k');
legend('5min','predicted','lowpass','Hourly avg','residual','Hourly residual',' ','location','best')
plot([t2(1) t2(end)],[0 0],'color',[.5 .5 .5],'linewidth',0.5)
ylabel('Sea level (m CD)')
title('DOT Port Geographe tide gauge')
xlabel([datestr(mean(Htime),'YYYY') ' (UTC+8)'])
grid on
box on
ax=gca;
ax.YTick=[-.5:.1:2.5];
ax.XTickLabelRotation=90;
ax.XLim=[t2(1) t2(end)];
ax.XTick=[floor(t2(1)):1:ceil(t2(end))];
datetick('x','mm/dd','keeplimits','keepticks')
ax.FontSize=14;

% ax.XLim=[datenum(2022,5,20) datenum(2022,5,27)]; % may storm zoom

% fname=['/Users/00068592/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/20220818_RBR_ABB_DUNS_EB_StuBarr/DOT_PortGeographe_tide_gauge_' datestr(t2(1),'yyyymmdd') '-' datestr(t2(end),'yyyymmdd')];
% fname=['/Users/00068592/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/20220818_RBR_ABB_DUNS_EB_StuBarr/DOT_PortGeographe_tide_gauge_' datestr(ax.XLim(1),'yyyymmdd') '-' datestr(ax.XLim(end),'yyyymmdd')];
fname=[whereput 'DOT_PortGeographe_tide_gauge_' datestr(t2(1),'yyyymmdd') '-' datestr(t2(end),'yyyymmdd')];

savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the spectral analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; % this is a dud figure created so it doesn't appear by default in previous figure;
% data=p(goodi); % for depth data (observed)
data=z2;

% data=zeros(2431001,1);
% data=zeros(length(P),1);

% data=AnDepthmm/1000.; % for conversion frommm to meters depth

% sample frequency 5 min = 300 sec
% sampleperiod=5*60;
% sampleperiod=3600;
% samplehertz=sample_rate_hz;
sampleperiod=300;

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

close(gcf)% this closes teh default figure that we don't want
disp('Spectra calculated')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To format the plot - Sxx vs f (log-log)

fh=figure;
fh.Color='w';
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
text(f24,ytop,' 24h','fontsize',[12])

% to add a vertical line at a period of 12 hours (f12)

f12=1/(12*60*60);
plot([f12 f12],ylim,'k:')
text(f12,ymiddle+260,' 12h','fontsize',[12])

% to add a vertical line at a period of 4 hours (f4)
f4=1/(4*60*60);
plot([f4 f4],ylim,'k:')
text(f4,ymiddle+100,' 4h','fontsize',[12])

% % 2hr
% f2=1/(1*60*60);
% plot([f2 f2],ylim,'k:')
% text(f2,ymiddle,' 1h','fontsize',[12])

% 15 mins
f3=1/((15/60)*60*60);
plot([f3 f3],ylim,'k:')
text(f3,ytop,'15 min','fontsize',[12])

% f3=1/((1/60)*60*60);
% plot([f3 f3],ylim,'k:')
% text(f3,ytop,'1 min','fontsize',[12])
% 
% f20=1/300;
% plot([f20 f20],ylim,'k:')
% text(f20,ymiddle,'300 sec','fontsize',[12])
% %
% f20=1/30;
% plot([f20 f20],ylim,'k:')
% text(f20,ymiddle,'30 sec','fontsize',[12])
% %
% f12=1/15;
% plot([f12 f12],ylim,'k:')
% text(f12,.2,'15 sec','fontsize',[12])
% title([instrument ' ' location],'interpreter','none')
% 
% f5=1/5;
% plot([f5 f5],ylim,'k:')
% text(f5,.6,'5 sec','fontsize',[12])
% title([instrument ' ' location],'interpreter','none')

% fname=['/Users/00068592/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/20220818_RBR_ABB_DUNS_EB_StuBarr/DOT_PortGeographe_tide_gauge_SPECTRA_' datestr(t2(1),'yyyymmdd') '-' datestr(t2(end),'yyyymmdd')];
fname=[whereput 'DOT_PortGeographe_tide_gauge_SPECTRA_' datestr(t2(1),'yyyymmdd') '-' datestr(t2(end),'yyyymmdd')];
savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TIME FREQ
% 
%   
% %         [logf logj timef]=get_tidespecft(t2,z2,4096,1024,1/300);
%         [logf logj timef]=get_tidespecft(t2,z2,16,4,1/300);
%         
%         % plot it
%         
%         
%         figure;
%         set(gcf,'Position',[200 100 900 500])
%         set(gcf,'PaperPositionMode','auto')
%         pcolor(timef,logf,logj(:,1:length(timef)))
%         shading('interp')
%         
%         hold on
%         
%         %axis([735600 735965  -12 -6])
%         % to add a vertical line at a period of 24 hours (f24)
%         ymiddle=0.20;
%         xlim=get(gca,'xlim');
%         
%         f24=log(1/(12*60*60));
%         plot(xlim,[f24 f24],'k-')
%         text(timef(100),f24+0.25,' 12h','fontsize',[12])
%         
%         f24=log(1/(24*60*60));
%         plot(xlim,[f24 f24],'k-')
%         text(timef(100),f24+0.25,' 24h','fontsize',[12])
%         
%         % f24=log(1/(6*60*60));
%         % %plot(xlim,[f24 f24],'k-')
%         % text(timef(100),f24+0.25,'6hr','fontsize',[12])
%         
%         % f24=log(1/2);
%         % plot(xlim,[f24 f24],'k-')
%         % text(timef(100),f24+0.25,' 2s','fontsize',[12])
%         
%         % f24=log(1/4.5);
%         % plot(xlim,[f24 f24],'k-')
%         % text(timef(100),f24+0.25,' 4.5s','fontsize',[12])
%         
%         f24=log(1/8);
%         plot(xlim,[f24 f24],'k-')
%         text(timef(100),f24+0.25,' 8s','fontsize',[12])
%         
%         f24=log(1/15);
%         plot(xlim,[f24 f24],'k-')
%         text(timef(100),f24+0.25,' 15s','fontsize',[12])
%         
%         f24=log(1/30);
%         plot(xlim,[f24 f24],'k-')
%         text(timef(100),f24+0.25,' 30s','fontsize',[12])
%         
%         % f24=log(1/130);
%         % plot(xlim,[f24 f24],'k-')
%         % text(timef(1000),f24+0.25,'130s','fontsize',[12])
%         
%         % f24=log(1/200);
%         % plot(xlim,[f24 f24],'k-')
%         % text(timef(1000),f24+0.25,'200s','fontsize',[12])
%         
%         
%         f24=log(1/300);
%         plot(xlim,[f24 f24],'k-')
%         text(timef(100),f24+0.25,' 300s','fontsize',[12])
%         %
%         % f24=log(1/600);
%         % plot(xlim,[f24 f24],'k-')
%         % text(timef(100),f24+0.25,' 600s','fontsize',[12])
%         
%         xlabel('Time')
%         ylabel('Frequency(Hz)')
%         
%         colormap(jet)
%         title([instrument ' ' location],'interpreter','none')
%         set(gca,'xlim',[timef(1) timef(end)],'xtick',xt,'fontsize',font_size)
%         datetick('x','mm/dd','keeplimits','keepticks')
%         rotateXLabels( gca(), 90 )
%         % set(gca,'xlim',[timef(1) timef(end)],'fontsize',font_size);
%         
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               save <1 hr sealevel data to text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% table
timestr=datestr(mtime,'yyyy-mm-dd HH:MM:SS');
Time_WST=cellstr(timestr);

    T = table(timestr,z,tide,residual);

outf_tname=[whereput 'DOT_PortGeographe_tide_gauge_'  datestr(mtime(1),'yyyymmdd') '-' datestr(mtime(end),'yyyymmdd') '_5_minutes.csv'];
writetable(T,outf_tname,'Delimiter',',','QuoteStrings',true)

disp('5min sea level data saved to csv text file')