


% 
% opts = detectImportOptions(filename)
% 
% getvaropts(opts,'Timestamp_AWST_')
% 
% 
% filename='Stirling_test.csv'

filename='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2020_Stirling/StirlingData2020.csv'

t=readtable(filename);
ts=datenum(t.Timestamp_AWST_);
datestr(ts,'yyyy-mm-dd HH:MM:SS.FFF');

figure; 
plot(ts+12/24,t.UnderwaterPressure_m_H20_-10,'color',[.7 .7 .7]);
datetick;
hold on

ax1=gca;


fm=load('/Users/00068592/Drop/YASHA (3)/1_DATA/RECENT/WADOT/fm.mat'); % UTC time
plot(fm.ts+8/24,fm.wl+2.8,'k')
ax2=gca;

load fmrawtg2020  %local time % time Sealevel


plot(time,Sealevel+2.4,'k','linewidth',2)
ax3=gca;


set(gca,'ylim',
[-2 2];



