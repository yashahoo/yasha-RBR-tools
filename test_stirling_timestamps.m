addpath(genpath(pwd))

filename='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/20210101_Stirling_test.txt'
Stirlingtest = importfile6(filename);

% Stirlingtest.time1.Format = 'yyyy-MM-dd';
% Stirlingtest.time3.Format = 'yyyy-MM-dd HH:mm:ss.SS';
t2=datenum(Stirlingtest.time1)+datenum(Stirlingtest.time3)-datenum(datestr(datenum(now),'yyyy-mm-dd'));

% datestr(t2(1:10),'yyyy-mm-dd HH:MM:SS.fff')
press=Stirlingtest.press;
press=press-mean(press,'omitnan');
press(press<-2)=NaN;



fm=load('/Users/00068592/Drop/YASHA (3)/2_ANALYSIS/3_MAT_FILES/fm.mat');

id=find(fm.ts>datenum(2020,1,1));
zt=fm.ts(id);
zz=fm.wl(id);

% Tpred=[datenum(2020,1,1):datenum(0,0,0,0,5,0):datenum(2021,2,1)];
Tpred=[datenum(2021,10,22):datenum(0,0,0,0,5,0):datenum(2021,11,1)];

[names,freq,tidecon,Htide]=t_tide(zz,'interval',1,'start',zt(1),'synthesis',1);
% mYout=t_predic(mts_sub,mnames,mfreq,mtidecon,'synthesis',1);
tide=t_predic(Tpred,names,freq,tidecon,'synthesis',1);

figure; plot(Tpred,tide);
datetick
hold on; plot(zt,zz);

 

% plot(t2,press);


t3=t2-20/24;
plot(t3,press);


t2=(datenum(Stirlingtest.time1)-1)+((datenum(Stirlingtest.time3))-datenum(datestr(datenum(now),'yyyy-mm-dd')));

figure; plot(Tpred,tide);
datetick
hold on; plot(zt,zz);

plot(t2,press);