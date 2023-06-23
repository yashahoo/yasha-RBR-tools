function WindArrows6(u,v,wtime,begintime,numdays,interval,xtickhrs,ymax,ymax2)


%plot seabreeze-like plot of wind arrows
% 
% yasha hetzel 2016
% 
%
% usage:
% 
% WindArrows6(u,v,wtime,begintime,numdays,interval,xtickhrs,ymax,ymax2)
% inputs:
% 
% -wind speed east/north components u,v in m/s
% -wind time in matlab datenum format
% -begintime = matlab string e.g. 'yyyy-mm-dd HH:MM'
% -numdays = number of days to plot, 
% -interval = spacing between arrows (in hours)
% (interval of 1 recommended for up to 10 days, 3 recommended for 1 month time period; looks bad longer than this)
% -xtickhrs = xtick spacing (in hours)
% -ymax = y axis upper limit,

% NOTE: color of arrows is determined by jet colormap with max color equal to ymax
%       be careful of setting numdays compared to start time to get desired
%       end date (may need to set to start of following day)
% required: arrows.m, cart2polCOMPASS.m
% 
% example:
% WindArrows4(u,v,wtime,'2009-6-1 00:00',7,1,12,20)
% this interpolates to hourly data using interp1, plots for 7 days, with
% ylim of 20 m/s, and 12 hour spacing between tick marks
%
% after loading  saved bom station data:
% si=20; when=datenum(2018,7,1); WindArrows4(archive(si).u,archive(si).v,archive(si).mtime_UTC,when,30,3,5*24,15); title(archive(si).name(1)); set(gcf,'color','w')
%
% revisions:
% added ability to choose interval between arrows yh 27/04/16
% updated color to map to 'jet' colormap yh 22/01/17
% added ability to set xtick spacing 'xtickhrs' yh 23/01/17
% updated help section yh 18/12/2018
% does not open new figure automatically so can be used in other axes yh 2019-05-23


% %required
% numdays=7;
% timestart=datenum(2009,6,1,0,0,0);
% %optional
% ymax=20;

SHAPE=[.6,.5,0.5,.25]; % SHAPE = [HEADW,HEADL,HEADI,LINEW] 


timestart=datenum(begintime);

% timefinish=datenum(2009,6,8,0,0,0);
timefinish=timestart+numdays;




% interpolate to nice 3 hourly time
ti=[timestart:interval/24:timefinish]; %datenum(time(1):1/24:time(end-1));
u=interp1(wtime,u,ti);
v=interp1(wtime,v,ti);
% rename
wtime=ti;



% better to use u,v as inputs to avoid confusion - convert to speed and
% direction
[dir spd]=cart2polCOMPASS(u,v);

% find limits
ti_start =find(wtime> timestart-2/24,1,'first');
ti_end=find(wtime< timefinish+2/24,1,'last');

% some stufff for plotting setup
dt=mean(diff(wtime))*24; % time between observations
dw=xtickhrs/dt; %24/dt; % 12 hr ticks vs 24 hr ticks
w=[1:length(wtime)];
wtick=w(1:dw:end)';
%wticklab=datestr(wtime(1:dw:end),'mm/dd HH:MM');
wticklab=datestr(wtime(1:dw:end),'ddd dd-mmm ');




% Black and white version
% figure;
% set(gcf,'Position',[200 100 1500 400])
% set(gcf,'PaperPositionMode','auto')
% hold on
% % plot(1:length(wtime),wspd/max(wspd),'k-')
% plot(w,spd,'k-')
% arrows(w,spd,2,dir,SHAPE,'FaceColor','w','EdgeColor','k','LineWidth',0.01)
% axis equal
% set(gca,'xlim',[ti_start ti_end],'XTick',wtick,'XTickLabel',wticklab,'xgrid', 'on') %,'XTickLabel',[],'Ytick',[] [w(1) w(end)] %'xlim',[1 numdays*24]
% set(gca,'ylim',[0 ymax])
% box on
% grid on
% xlabel(['Date ' '(' datestr(mean(wtime),'yyyy') ')'])
% ylabel('Wind speed (m s^-^1)')



%% plot loop with color
% [option to plot with color]

 cm=colormap(jet);
 
% figure;
% set(gcf,'Position',[200 100 1500 400])
% set(gcf,'PaperPositionMode','auto')
hold on
plot(w,spd,'k-')

for i=1:length(wtime)
    
    sp=spd(i);
    t=w(i);
    d=dir(i);

%     get the right color
    rat=sp/ymax;
    try
    color=cm(round((rat).*length(cm)),:);
    catch
        
        if rat==0
        color=cm(1,:);
        
        elseif rat>1
            color=cm(length(cm),:);
            
         elseif isnan(rat)
             color=[1 1 1]; % make white if nan
            
        end
    end
    arrows(t,sp,2,d,SHAPE,'FaceColor',color,'EdgeColor','k','LineWidth',0.01)
%     arrows(t,sp,2,d,SHAPE,'FaceColor',[sp/max(spd) 0   1-sp/max(spd) ],'EdgeColor','k','LineWidth',0.01)
    axis equal
    
end

%set(gca,'xlim',[ti_start ti_end+3],'XTick',wtick,'XTickLabel',wticklab,'xgrid', 'on') %,'XTickLabel',[],'Ytick',[] [w(1) w(end)]
set(gca,'xlim',[ti_start ti_end],'XTick',wtick,'XTickLabel',wticklab,'xgrid', 'on')
if exist('ymax2','var')
set(gca,'ylim',[0 ymax2])
else
set(gca,'ylim',[0 ymax])
end

grid on
box on
xlabel(['Date ' '(' datestr(mean(wtime),'yyyy') ')'])
%ylabel('Wind speed (m s^-^1)')
%ylabel('Wind speed (knots)')
%set(gca,'fontsize',14)
%% junk
% clc;
% clear all;
% close all;
% 
% load Butterfilt_winds_Ocenography_2009_2010_2011-v2.txt
% 
% %--loading wind data
% WTIME=Butterfilt_winds_Ocenography_2009_2010_2011_v2(:,1); 
% WinSpd=Butterfilt_winds_Ocenography_2009_2010_2011_v2(:,4);
% WinDir=Butterfilt_winds_Ocenography_2009_2010_2011_v2(:,5);

% wspd=interp1(wtime,wspd,ti);
% wdir=interp1(wtime,wdir,ti);

% reaname to fit code
% wspd=WinSpd;
% wdir=WinDir;
% wtime=WTIME;
% [u v]=pol2cartCOMPASS(wdir,wspd);

% better to use u,v as inputs to avoid confusion
% [u v]=pol2cartCOMPASS(wdir,wspd);
% [dir spd]=cart2polCOMPASS(-u,-v);



%old way
% % subset data (hourly preferred)
% % find avg tiem between obs and find how many to skip??
% ddt=round(1/(mean(diff(WTIME))*24));
% % keep only hourly
% wspd=WinSpd(1:ddt:end);
% wdir=WinDir(1:ddt:end);
% wtime=WTIME(1:ddt:end);





