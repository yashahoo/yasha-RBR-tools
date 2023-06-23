% yasha hetzel April 2021
% this m-file reads raw .rsk file downloaded from RBR pressure sensors
% using rsk-tools (downloadad from rbr website)
% and saves to netcdf file. then it reads the netcdf file and calculates
% waves (optional), tidal analysis, lowpass filter, highpass filter (optional) 
% and makes a number of plots that are saved to a directory along
% with the processed data as matfile (optional) or the hourly wave data as a csv text file
% **There are a lot of m-files required so you need to set the path to make 
% sure all the m-files in the yasha_RBR_tools are available.
%%


clear all;
clc;
close all;

% required:
% t-tide
% rotateXlabels
% subaxis
% exportfig
% lanczosfilter
% and some others including rbr rsk-tools


restoredefaultpath % this just fixes if there are any confliting functions

% set the paths here relevant to your computer
% addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/yashatools/yasha_RBR_tools')); % [edit]  this should contain all the functions you need, except the rsk-tooks that are required to read the raw data (see below)
% % addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rbr-rsktools-2.3.0')); %  [edit] this should be installed on your computer somewhere (it has been tested ususing old version)
% addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rbr-rsktools')); % use this old version to ensure it reads temperature!

% this should work if you run it from the yasha_RBR_tools directory
addpath(genpath(pwd))
addpath(genpath('rbr-rsktools'))
%%-------------------------------------------------------------------------%
%%  ------------------------settings-------------------------------------- %
%%-------------------------------------------------------------------------%
% set location and name of raw file. by default, figures etc will be saved
% to same directory
% location='EastBusselton'; % used for filenaming
% rawfile='/Users/00068592/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/SN206857_Dirk_Port_Geographe_20210831-20211021/206857_20220629_0954.rsk'
% location='Abbey'; rawfile='/Users/00068592/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/20220818_RBR_ABB_DUNS_EB_StuBarr/ABBEY_206860_20220818_0948.rsk'
% location='OldDunsborough'; rawfile='/Users/00068592/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/20220818_RBR_ABB_DUNS_EB_StuBarr/DUNS_206857_20220818_1312.rsk'
location='EastBusselton'; rawfile='/Users/00068592/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/20220818_RBR_ABB_DUNS_EB_StuBarr/EASTBUSSO_209847_20220818_1102.rsk'

% deployment info
sample_rate_hz=2; % sampling in hz (normally 2)
max_period=25; % seconds for cutoff in calculating waves
subsample_mins=1; % water level averaging interval in minutes (all other calculations done in hourly intervals)

% exta options -   'Y' if you want to calculate/plot; 'N' if don't
round_start_finish = 'Y'       % if you want nice start and end times (e.g. on the 1/2 hour) 
calculate_waves    = 'Y'       % calculate wave parameters including Hs, Tz, etc.
calculate_IG       = 'Y'       % calculate Infragravity waves (shoudl be 'Y' if calculating waves)
get_time_freq      = 'Y'       % use chari's function to calculate time-frequency plot
do_highpass_filter = 'N'       % do highpass filter to remove waves > max_period
plot_highpass_weeks= 'N'       % plot weekly highpass filter (was useful to look at boatwakes in albany)
save_matfile       = 'N'       % save the processed data to matfile
save_matfig        = 'N';      % save the .fig files?
print_fig          ='Y';       % print to image file?
use_export_fig     ='N'        % if print_fig ='Y' -->  'N' uses native matlab print function (default) OR 'Y' uses export_fig to make nice plots
img_type           = {'pdf','png'};  % if print_fig ='Y' -->  {'pdf','png','eps'} % options. can be multiple. pdf will not save if native matlab unless edit function 


% set where to save figures and data
fs=12; % set font size for figs
whereput=[rawfile(1:end-4) '_processed'];
mkdir(whereput);


% time range to process. If unknown, make NaN and will be selected interactively
% begin=datenum(2020,7,23,16,23,0);
% finish=datenum(2020,9,27,18,0,0);
begin=NaN;
finish=NaN;

% flag bad data with nans [if no bad data make badstart/end = NaN]. this
% seems to cause some problemss,so avoid for now if possible
% badstart=datenum(2018,2,4,4,17,3);
% badend=datenum(2018,2,4,4,19,20);
badstart=nan;
badend=nan;

%%-------------------------------------------------------------------------%
%%-------------------------------------------------------------------------%
%%-------------------------------------------------------------------------%
%% read and create the nc file and begin processing
tic
if isempty( regexp( rawfile, filesep, 'start' )) ;
    nopath=1;
else
    nopath=0;
end

if nopath>0
    filebase=rawfile(1:end-4);
else
naminds=regexp( rawfile, filesep, 'start' ) ;
end
filebase=rawfile(naminds(length(naminds))+1:length(rawfile)-4);
ncref=[whereput filesep filebase '.nc'];



% check that nc file doesnt yet exist, delete if it does
if isfile(ncref) %exist(ncref,'file') %exist('ncref','var')
    disp('Netcdf output file already exists, delete it?')

    promptMessage = sprintf('Overwrite netcdf output file and CONTINUE or QUIT?');
    promptTitle = sprintf('Netcdf file exists');
    button = questdlg(promptMessage,promptTitle,'Quit','Continue','Continue');
    if strcmpi(button, 'Quit')
        close all
        return; % Or break or continue
    elseif strcmpi(button, 'Continue')
        disp(['Deleting: ' ncref])    
        delete(ncref)
        clear ncref  
    end  
%     disp(['Deleting: ' ncref])
    
end

% read raw file and save to netcdf
[ncref]=read_RBRyh(rawfile,location,whereput);
% ncref='/Volumes/Margs_Clone/RBRpressure_sensors/SN202078_Smiths_20210409/202078_20210409_1525_processed/202078_20210409_1525.nc'

%%
% % % this used later
% ncref=[whereput '/' rawfile(end-23:end-4) '.nc'];
close all;


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

% if you want nice start and end times (e.g. on the 1/2 hour) 
switch round_start_finish
    case 'Y'
        tv_dv_begin = datevec(begin);tv_dv_finish = datevec(finish);
        begin = datenum(datevec(datenum([tv_dv_begin(:,1:4) [30*(tv_dv_begin(:,5)<30)+60*(tv_dv_begin(:,5)>=30)] 0])));  % <30 --> 30
        finish = datenum(datevec(datenum([tv_dv_finish(:,1:4) [0*(tv_dv_finish(:,5)<30)+30*(tv_dv_finish(:,5)>=30)] 0])));  %  rounds down to nearest half hour
end

% find data before/after deployment / retrieval
badi=find(ts<begin | ts>finish); %

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
clear Hs Ts Htime Hrms Hmax Tmax H10 T10 Tz psub tsub inds


sph=sample_rate_hz*3600;

switch calculate_waves
    case 'Y'
        disp('Calculating HOURLY waves + sealevels')
    case 'N'
        disp('NO waves calculated')
        disp('Calculating HOURLY sealevels')
end

numdays=floor(ts(end)-ts(1)); % yh ?? -1
clear waves
co=0;
for hr=1:numdays*24
   co=co+1; 
   inds=[hr*sph-(sph-1):hr*sph];
   if max(inds)>length(p) % yh 20210128 test      
       disp('NOT CALCULATING LAST HOUR!!')
       break
   end
psub=p(inds);
tsub=ts(inds);

switch calculate_waves
    case 'Y'
        clear waves
        waves=wavepar_yh(psub,sample_rate_hz,1/max_period); %2 hz
        
        Hs(co)=waves.Hs;
        Ts(co)=waves.Ts;
        Hrms(co)=waves.Hrms;
        Hmax(co)=waves.Hmax;
        Tmax(co)=waves.Tmax;
        H10(co)=waves.H10;
        T10(co)=waves.T10;
        Tz(co)=waves.Tz;

end

Htime(co)=tsub(sph/2+1);
Hpress(co)=nanmean(psub); % hourly pressure

if exist('temp','var')
    tempsub=nanmean(temp(inds));
Htemperature(co)=tempsub;
end

end

if exist('Htemperature','var')
% if incomplete waves cut the temperature to be the same
if length(Htemperature)>length(Htime)
    Htemperature(end)=[];
    disp('removed last hour of Temperature data!!')
end
end
xlims=[ts(1) ts(end)]; % all
xlims=[begin finish];% good
xt=[xlims(1):xlims(2)];
yt=[0:2:30];




switch calculate_waves
    case 'Y'
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
        set(gca,'ylim',[0 ceil(max(Hmax))])
        datetick('x','mm/dd','keeplimits','keepticks')
        rotateXLabels( gca(), 90 )
        xlabel(datestr(mean(Htime),'YYYY'))
        title(instrument)
        fname=[whereput '/' instrument '_WAVE_HEIGHT_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
        
        figure;
        set(gcf,'Position',[200 100 700 300])
        set(gcf,'PaperPositionMode','auto')
        plot(Htime,T10)
        hold on
        plot(Htime,Ts)
        set(gca,'xlim',xlims,'xtick',xt,'ytick',yt) 
        plot(Htime,Tz)
        legend('T10','Ts','Tz')
        ylabel('Period (s)')
        datetick('x','mm/dd','keeplimits','keepticks')
        xlabel(datestr(mean(Htime),'YYYY'))
        rotateXLabels( gca(), 90 )
        fname=[whereput '/' instrument '_WAVE_PERIOD_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
        
end
 
 if exist('temp','var')
     figure;
     set(gcf,'Position',[200 100 700 300])
     set(gcf,'PaperPositionMode','auto')
     plot(Htime,Htemperature)
     hold on
%      legend('Hs','H10','Hmax')
     ylabel('temperature ( deg C)')
     set(gca,'xlim',xlims,'xtick',xt)
     set(gca,'ylim',[floor(min(Htemperature)) ceil(max(Htemperature))])
     datetick('x','mm/dd','keeplimits','keepticks')
     xlabel(datestr(mean(Htime),'YYYY'))
     rotateXLabels( gca(), 90 )
     title(instrument)
     fname=[whereput '/' instrument '_TEMPERATURE_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
 end
 
%% subsample to less than 1 hr



numdays=floor(ts(end)-ts(1)); % yh ?? -1
clear inds spw tsub psub time_subsampled press_subsampled
co=0;
for sub_hr=1:numdays*24*(60/subsample_mins)
   co=co+1; 
   spw=sph.*subsample_mins/60 ;  % samples per window sph=sample_rate_hz*3600;
   inds=[sub_hr*spw-(spw-1):sub_hr*spw];
   if max(inds)>length(p) % yh 20210128 test      
       disp('NOT CALCULATING LAST time slot!!')
       break
   end
psub=p(inds);
tsub=ts(inds);


time_subsampled(co)=tsub(spw/2+1);
press_subsampled(co)=nanmean(psub); % hourly pressure

if exist('temp','var')
    tempsub=nanmean(temp(inds));
temperature_subsampled(co)=tempsub;
end

end

 figure; 
 plot(time_subsampled,press_subsampled)
datetick('keeplimits')
title(['subsample interval = ' num2str(subsample_mins) ' minutes']); 
 
 %% get IG waves (this overlaps with above)

switch calculate_waves
    case 'Y'
        
        [wavedata hs1  hm1 S]=get_chari_wave_analysis(ts,p,2); %[wavedata hs1  hm1 S]=get_chari_wave_analysis(time,data,srate); srate in Hz 2/second
        
        figure;
        set(gcf,'Position',[200 100 700 300])
        set(gcf,'PaperPositionMode','auto')
        hold on
        box on
        plot(wavedata(:,1),wavedata(:,25))
        set(gca,'xlim',xlims,'xtick',xt)
        ylabel('IG height (m)')
        datetick('x','mm/dd','keeplimits','keepticks')
        rotateXLabels( gca(), 90 )
        xlabel(datestr(mean(Htime),'YYYY'))
        % figure;
        %  plot(wavestata(1,:),wavestata(25,:))
        %  hold on
        %  xlabel('Time (days)')
        %  ylabel('IG height (m)')
        fname=[whereput '/' instrument '_IG_WAVES_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
        
        
        
        % get hourly ig
        tig=wavedata(:,1);
        ig=wavedata(:,25);
        
        Hig=interp1(tig, ig,Htime);
        
end
%% get lowpass sea level (lanczos or butterworth 38 hours)

% SN124028
% p(1:289286)=0;
% p(9325510:end)=0;


% if incomplete waves cut the temperature to be the same
if length(Hpress)>length(Htime)
    Hpress(end)=[];
    disp('removed last hour of Hpress data!!')
end

% pad with zeros to make sure start is ok
Hpress0=Hpress;
% Hpress0(1:8)=0;
% Hpress0(1294:end)=0;
Hpress0=[zeros(1,100) Hpress0 zeros(1,100)]; 

low_cutoff=38;
% lowpass filter model (different method- not used)
[lp_wl,B,A]=lp_filter(6,low_cutoff,1,Hpress0);
Hpress0(1:100)=[];Hpress0(end-99:end)=[];
lp_wl(1:100)=[];lp_wl(end-99:end)=[];
lp_wl(1:28)=nan;
lp_wl(1029:end)=nan;
Hpress0(1:28)=NaN;
Hpress0(1029:end)=NaN;
% 
% figure; 
% hold on;
% plot(Htime,Hpress0)
% plot(Htime,lp_wl)
% datetick('x','yyyy/mm/dd HH:MM:SS','keeplimits')
% 
% figure; 
% hold on;
% plot(Htime,Hpress-lp_wl)
% % plot(Htime,lp_wl)
% datetick('x','yyyy/mm/dd HH:MM:SS','keeplimits')
% 



% try lanczos filter
      t=Htime;
      xn=Hpress;
      dT=1;
      Tc = 38; % hr
  [lanczos_wl,c,h,Cx,f] = lanczosfilter(xn,dT,1/Tc,[],'low');

%   lanczos_wl=xs;
%   
  
% choose which to use (they are the same from both filters)
lowpass=lanczos_wl-nanmean(lanczos_wl);

disp([num2str(Tc) ' hour Lanczos lowpass filter applied'])
% lowpass(badi)=nan;
% lowpass=lowpass-nanmean(lowpass);


% badi=[];
% badi=find(Htime<datenum(2017,9,22,0,0,0) | Htime>datenum(2017,11,11,11,0,0)); %SN124096

%% get tidal water levels by subtracting lowpass

% Htide=Hpress-lp_wl;

%% or do t-tide to get tides

     % predict  tide & get residual from model data
        clear names freq tidecon Htide

[names,freq,tidecon,Htide]=t_tide(Hpress,'interval',1,'start',Htime(1),'synthesis',1);
% mYout=t_predic(mts_sub,mnames,mfreq,mtidecon,'synthesis',1);

% get higher frequency tide output matching subsampled
% time_subsampled press_subsampled
tide_subsampled=t_predic(time_subsampled,names,freq,tidecon,'synthesis',1);
% calculate non-tidal residual
residual_subsampled=press_subsampled-nanmean(press_subsampled)-tide_subsampled;

f=figure; 
f.Position=[200 100 900 400];
hold on
plot(time_subsampled,press_subsampled-nanmean(press_subsampled),'k')
plot(time_subsampled,tide_subsampled,'b')
plot(time_subsampled,residual_subsampled,'r')
legend('total','tide','residual','location','best')
ylabel('Water level (m)')
xlabel(['Date (' datestr(time_subsampled(1),'YYYY') ')'])
set(gca,'xlim',xlims,'xtick',xt)
datetick('x','mm/dd','keeplimits','keepticks')
rotateXLabels( gca(), 90 )
box on;
grid on;
title(['Averaging interval = ' num2str(subsample_mins) ' minutes']); 

fname=[whereput '/' instrument '_SEALEVEL_' num2str(subsample_mins) 'mins_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
        
        
        
% calculate non-tidal residual
residual=Hpress-nanmean(Hpress)-Htide;

% PLOT all sea levels
f=figure; 
f.Color='w';
set(gcf,'Position',[200 100 900 400]) 
set(gcf,'PaperPositionMode','auto')
hold on
plot(Htime,Hpress-nanmean(Hpress),'k');
plot(Htime,Htide)
plot(Htime,residual,'r')
% plot(Htime,lp_wl-nanmean(lp_wl),'g')
plot(Htime,lanczos_wl-nanmean(lanczos_wl),'g');
% legend('Total  water level','Tide')
ylabel('Sea level (m)')
xlabel(datestr(mean(Htime),'YYYY'))
legend('Total','Tide','Residual','lowpass','location','best') % ,'lanczos lowpass'
set(gca,'xlim',xlims,'xtick',xt)
datetick('x','mm/dd','keeplimits','keepticks')
rotateXLabels( gca(), 90 )
box on;
grid on;
fname=[whereput '/' instrument '_SEALEVEL_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
 
 disp('T_tide analysis completed')

%% plot all

switch calculate_waves
    case 'Y'
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
            
            f=figure;
            f.Position=[200 100 600 800];
            f.Color ='w';
            subaxis(5,1,1)
            plot(Htime,Hmax)
            hold on;
            plot(Htime,Hs)
            plot(Htime,Hrms)
            legend('Hmax','Hs','Hrms')
            set(gca,'xlim',xlims,'xtick',xt)
            set(gca,'ylim',[0 ceil(max(Hmax))],'ytick',[0:.5:5])
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
            xlabel(datestr(mean(Htime),'YYYY'))
            
            
            
            
            fname=[whereput '/' instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
            
        end
        
        case 'N'
        disp('NO waves plotted')
end

%% save Hourly data to text

switch calculate_waves
    case 'Y'
        
        Hs=Hs';Sealevel=Hpress'; Hrms=Hrms';H10=H10';Hmax=Hmax'; Hlowpass_wl=lowpass; %lp_wl';
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
        outf_tname=[whereput '/' instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '_waves_hourly.csv'];
        writetable(T,outf_tname,'Delimiter',',','QuoteStrings',true)
        

        disp('Hourly wave and sea level data saved to csv text file')
        
end


switch calculate_waves
    case 'N'
        
        Sealevel=Hpress'; Hlowpass_wl=lowpass; %lp_wl';
        Htide=Htide'; Hresidual=residual';
         
        
        if exist('Htemperature','var')
            Htemperature=Htemperature';
        end
        % varnames={'Time WST','Hrms',
        
        % table
        timestr=datestr(Htime','yyyy-mm-dd HH:MM:SS');
        Time_WST=cellstr(timestr);
        if exist('Htemperature','var')
            T = table(timestr,Sealevel,Hlowpass_wl,Htide,Hresidual,Htemperature);
        else
            T = table(timestr,Sealevel,Hlowpass_wl,Htide,Hresidual);
        end
        % % outf_tname=[whereput '/' instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '.csv'];
        % outf_tname=[dir_ref instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '.csv'];
        outf_tname=[whereput '/' instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '_Sealevel_no_waves_hourly.csv'];
        writetable(T,outf_tname,'Delimiter',',','QuoteStrings',true)
        
         disp('Hourly sea level only data saved to csv text file')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save <1 hr sealevel data to text file




sub_time=time_subsampled';
sub_tide=tide_subsampled';
sub_residual=residual_subsampled';
sub_press=press_subsampled';

        if exist('temperature_subsampled','var')
            sub_temp=temperature_subsampled;  
        end
        % varnames={'Time WST','Hrms',
        
        % table
        timestr=datestr(sub_time,'yyyy-mm-dd HH:MM:SS');
        Time_WST=cellstr(timestr);
        if exist('sub_temp','var')
            T = table(timestr,sub_press,sub_tide,sub_residual,sub_temp);
        else
            T = table(timestr,sub_press,sub_tide,sub_residual);
        end
        % % outf_tname=[whereput '/' instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '.csv'];
        % outf_tname=[dir_ref instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '.csv'];
        outf_tname=[whereput '/' instrument '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '_Sealevel_no_waves_' num2str(subsample_mins) '_minutes.csv'];
        writetable(T,outf_tname,'Delimiter',',','QuoteStrings',true)
        
         disp('sub-Hourly sea level only data saved to csv text file')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the spectral analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; % this is a dud figure created so it doesn't appear by default in previous figure;
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
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)


%% use chari's function to calculate time-frequency plot

switch get_time_freq
    case 'Y'
        
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
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
     
        case 'N'
        disp('NO time-frequency plot created')
        
end

% 'VerticalAlignment', 'top'

 %% %% do highpass filter
 
 switch do_highpass_filter
     case 'Y'
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
                     % print('-dpng','-r200', fname);
                     % print('-dpdf', '-painters', fname);
                     %         export_fig(fname,'-pdf','-png','-transparent')
                     savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
                     close
                     
                 end
         end


         
         figure;
         hold on
         % set(gcf,'Position',[200 100 1200 800])
         set(gcf,'Position',[200 100 1200 400])
         set(gcf,'PaperPositionMode','auto')
         plot(tsub,data);
         % set(gca,'xlim',[ts(1) ts(end)]);
         % set(gca,'xticklabel',[]);
         % datetick('x','yyyy/mm/dd ','keeplimits')
         plot(wt(tf),wz(tf)-1,'r.')
         set(gca,'xlim',xlims,'xtick',xt)
         set(gca,'ylim',[floor(min(data)) ceil(max(data))],'ytick',[-1:.25:2])
         datetick('x','mm/dd','keeplimits','keepticks')
         rotateXLabels( gca(), 90 )
         ylabel('Height (m)')
         title([instrument ' highpass (<' num2str(max_period) ' s) sea level'])
         set(gca,'fontsize',12)
         box on
         grid on
         
         fname=[whereput '/' instrument '_HIGHPASS_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
         % print('-dpng','-r200', fname);
         % print('-dpdf', '-painters', fname);
         %         export_fig(fname,'-pdf','-png','-transparent')
         savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
         
     case 'N'
         disp('NO Highpass filter applied')
         
end
%% save all data to mat file

switch save_matfile
    case 'Y'
% fname=[dir_ref instrument '_PROCESSED_DATA_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];

fname=[whereput '/' instrument '_PROCESSED_DATA_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
save(fname)
disp(['All processed variables saved to: ' fname '.mat'])

    case 'N'
        disp('NOT saving matfile')
end

%%

disp('DONE!!!!')
toc
