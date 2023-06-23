% compare Success data to FRE tide gauge and waves to get timezone
% yasha hetzel 2023-06-16


%% load fm tide gauge and predict tide


% load wind data;
% load bom winds
bom=load('BOM_obs_archived.mat');
b=bom.archive(19);
%%

% load tide gauge
fm=load('/Users/00068592/Drop/YASHA (3)/2_ANALYSIS/3_MAT_FILES/fm.mat');

id=find(fm.ts>datenum(2020,1,1));
tgt=fm.ts(id);
tgz=fm.wl(id);

% Tpred=[datenum(2020,1,1):datenum(0,0,0,0,5,0):datenum(2021,2,1)];
Tpred=[datenum(2021,10,22):datenum(0,0,0,0,5,0):datenum(2021,11,1)];

[names,freq,tidecon,Htide]=t_tide(tgz,'interval',1,'start',tgt(1),'synthesis',1);
% mYout=t_predic(mts_sub,mnames,mfreq,mtidecon,'synthesis',1);
tide=t_predic(Tpred,names,freq,tidecon,'synthesis',1);


%% do for 2021

% load the raw Success data
dat=load('/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2021-2023_Success_CLEANED_v1/202110_Cleaned_pressure_data.mat');
rt=dat.data(:,1);
rz=dat.data(:,2)-mean(dat.data(:,2),'omitnan');


pdat=load('/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2021-2023_Success_CLEANED_v1/202110_Cleaned_pressure_data_Success_processed_spectral/FPA_pressure_sensor_Success_PROCESSED_DATA_20211022-20211031.mat');
pt=pdat.sub_time;
pz=pdat.sub_press-mean(pdat.sub_press,'omitnan');

waves=pdat.Waves_output_Table;
wt=datenum(waves.timestr);
wz=waves.Hm0sea;

figure; 
hold on; 
plot(rt,rz);
plot(Tpred,tide,'k');
plot(pt,pz);
plot(wt,wz,'r');
datetick
legend('Success','FM tide gauge','5 min success','Hm0sea')
set(gca,'fontsize',14)
xlabel('Time (UTC)')
grid on


figure; 
hold on; 
plot(rt+8/24,rz);
plot(Tpred+8/24,tide,'k');
plot(pt+8/24,pz);
plot(wt+8/24,wz,'r');
plot(wt+8/24,waves.Hm0swell,'b')
datetick
legend('Success','FM tide gauge','5 min success','Hm0sea','Hm0swell')
set(gca,'fontsize',14)
xlabel('Time (AWST)')
grid on




% plot(b.mtime_UTC+8/24,b.wind_spd_kt./10)
plot(b.mtime_UTC+8/24,b.v/10)
legend('Success','FM tide gauge','5 min success','Hm0sea','Hm0swell','n-component of wind /10')
title('Success CombinedDataFeed')
box on

%% do for 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tide gauge prediction
Tpred=[datenum(2022,6,1):datenum(0,0,0,0,5,0):datenum(2022,7,1)];

% [names,freq,tidecon,Htide]=t_tide(tgz,'interval',1,'start',tgt(1),'synthesis',1);
% % mYout=t_predic(mts_sub,mnames,mfreq,mtidecon,'synthesis',1);
tide=t_predic(Tpred,names,freq,tidecon,'synthesis',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the raw Success data
dat=load('/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2021-2023_Success_CLEANED_v1/202206_Cleaned_pressure_data.mat');
rt=dat.data(:,1);
rz=dat.data(:,2)-mean(dat.data(:,2),'omitnan');

% processed wave data
pdat=load('/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2021-2023_Success_CLEANED_v1/202205_Cleaned_pressure_data_Success_processed_spectral/FPA_pressure_sensor_Success_PROCESSED_DATA_20211022-20211031.mat');
pt=pdat.sub_time;
pz=pdat.sub_press-mean(pdat.sub_press,'omitnan');

waves=pdat.Waves_output_Table;
wt=datenum(waves.timestr);
wz=waves.Hm0sea;

figure; 
hold on; 
plot(rt,rz);
plot(Tpred,tide,'k');
% plot(pt,pz);
% plot(wt,wz,'r');
datetick
legend('Success','FM tide gauge','5 min success','Hm0sea')
set(gca,'fontsize',14)
xlabel('Time (UTC)')
grid on


figure; 
hold on; 
plot(zt+8/24,zz);
plot(Tpred+8/24,tide,'k');
plot(pt+8/24,pz);
plot(wt+8/24,wz,'r');
plot(wt+8/24,waves.Hm0swell,'b')
datetick
legend('Success','FM tide gauge','5 min success','Hm0sea','Hm0swell')
set(gca,'fontsize',14)
xlabel('Time (AWST)')
grid on




% plot(b.mtime_UTC+8/24,b.wind_spd_kt./10)
plot(b.mtime_UTC+8/24,b.v/10)
legend('Success','FM tide gauge','5 min success','Hm0sea','Hm0swell','n-component of wind /10')
title('Success CombinedDataFeed')
box on


