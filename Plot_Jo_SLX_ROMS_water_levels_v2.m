% plot timeseries ecxtracted from schism and roms cwa for Jo 
% yasha hetzel June 2022


clear all; close all; clc;
% load SLX schism data
load('/Users/00068592/Documents/RESEARCH/MODELING/SCHISM/AUS/LONG_RUN/buckee/v8_COMBINED_sealevel_flux_figs_data/Jo_paper_sites_combined_data_v8.mat')

slx.pos=site_data(1).combined;

a=load('/Users/00068592/Documents/RESEARCH/MODELING/SCHISM/AUS/LONG_RUN/buckee/extracted_SCHISM_sealevel_components_11499.mat')
raw.pos.z=a.raw_model_sealevel;
raw.pos.time=a.time;
raw.pos.zeta=raw.pos.z./1000;
% raw.zeta=raw.z-mean(raw.z,'omitnan');


% load roms reanalysis data
ncref='/Users/00068592/GitHub/cwa-roms/POS_2001-2020_qck_zeta.nc'
roms.pos.time=ncread(ncref,'ocean_time')./3600/24 +datenum(2000,1,1);
roms.pos.zeta=squeeze(ncread(ncref,'zeta'));

ncref='/Users/00068592/GitHub/cwa-roms/POS_2020-202205_OPERATIONAL_qck_zeta.nc'
roms.pos.operational.time=ncread(ncref,'ocean_time')./3600/24 +datenum(2000,1,1);
roms.pos.operational.zeta=squeeze(ncread(ncref,'zeta'));

%load tide gauge data
gn=load('gn.mat');



%% save sealevels to csv files

whereput=[pwd filesep 'SealevelX_ROMS_sealevels_v2']
mkdir(whereput);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             slx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[years months days hours minutes seconds]=datevec(slx.pos.time);
site_id=slx.pos.site_id;
%  save to CSV comma delimted file
disp('Saving data to csv file...')
output_filename=[whereput '/' site_id '_sealevelX_sealevels']


    fid = fopen([output_filename '.csv'],'w');
%     
%     
    fprintf(fid,'Ivica SLX schism water levels for Jo Buckee study sites \n');
    fprintf(fid,'Created by Yasha Hetzel:               %s (UTC+8)\n',datestr(now));
    fprintf(fid,'With:                                  %s\n','plot_Jo_SLX_ROMS_water_levels.m');
    fprintf(fid,'Name:                                  %s\n','Post Office Island, Abrolhos WA');
    fprintf(fid,'ID:                                    %s\n',site_id);
    fprintf(fid,'Nearest grid point Longitude           %0.4f\n',slx.pos.lon);
    fprintf(fid,'Nearest grid point Latitude            %0.4f\n',slx.pos.lat);
    fprintf(fid,'Times are in UTC/GMT, To convert to local time add: %s hours\n',num2str(slx.pos.timezone));
    fprintf(fid,'\n');
    fprintf(fid,'Year, Month, Day, Hour, Total sea level (m)\n');
    fprintf(fid, '%d,%d,%d,%d,%0.4f\n',   [years, months, days, hours, slx.pos.total_adjusted_sealevel./1000]')        
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       roms reanalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[years months days hours minutes seconds]=datevec(roms.pos.time);
site_id=slx.pos.site_id;
%  save to CSV comma delimted file
disp('Saving data to csv file...')
output_filename=[whereput '/' site_id '_ROMS_REANALYSIS_sealevels']


    fid = fopen([output_filename '.csv'],'w');
%     
%     
    fprintf(fid,'Ivica Cwa-Roms REANALYSIS water levels for Jo Buckee study sites \n');
    fprintf(fid,'Created by Yasha Hetzel:               %s (UTC+8)\n',datestr(now));
    fprintf(fid,'With:                                  %s\n','plot_Jo_SLX_ROMS_water_levels.m');
    fprintf(fid,'Name:                                  %s\n','Post Office Island, Abrolhos WA');
    fprintf(fid,'ID:                                    %s\n',site_id);
    fprintf(fid,'Nearest grid point Longitude           %0.4f\n',slx.pos.lon);
    fprintf(fid,'Nearest grid point Latitude            %0.4f\n',slx.pos.lat);
    fprintf(fid,'Times are in UTC/GMT, To convert to local time add: %s hours\n',num2str(slx.pos.timezone));
    fprintf(fid,'\n');
    fprintf(fid,'Year, Month, Day, Hour, Total sea level (m)\n');
    fprintf(fid, '%d,%d,%d,%d,%0.4f\n',   [years, months, days, hours, roms.pos.zeta-mean(roms.pos.zeta,'omitnan')]')        
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       roms operational
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[years months days hours minutes seconds]=datevec(roms.pos.operational.time);
site_id=slx.pos.site_id;
%  save to CSV comma delimted file
disp('Saving data to csv file...')
output_filename=[whereput '/' site_id '_ROMS_OPERATIONAL_sealevels']


    fid = fopen([output_filename '.csv'],'w');
%     
%     
    fprintf(fid,'Ivica Cwa-Roms OPERATIONAL water levels for Jo Buckee study sites \n');
    fprintf(fid,'Created by Yasha Hetzel:               %s (UTC+8)\n',datestr(now));
    fprintf(fid,'With:                                  %s\n','plot_Jo_SLX_ROMS_water_levels.m');
    fprintf(fid,'Name:                                  %s\n','Post Office Island, Abrolhos WA');
    fprintf(fid,'ID:                                    %s\n',site_id);
    fprintf(fid,'Nearest grid point Longitude           %0.4f\n',slx.pos.lon);
    fprintf(fid,'Nearest grid point Latitude            %0.4f\n',slx.pos.lat);
    fprintf(fid,'Times are in UTC/GMT, To convert to local time add: %s hours\n',num2str(slx.pos.timezone));
    fprintf(fid,'\n');
    fprintf(fid,'Year, Month, Day, Hour, Total sea level (m)\n');
    fprintf(fid, '%d,%d,%d,%d,%0.4f\n',   [years, months, days, hours, roms.pos.operational.zeta-mean(roms.pos.operational.zeta,'omitnan')]')        
fclose(fid);

%%
% plot them
f=figure; 
f.Position=[69 439 2143 738];
f.Color='w';

% plot(raw.pos.time,raw.pos.zeta-mean(raw.pos.zeta,'omitnan'),'c');
hold on;
plot(slx.pos.time,slx.pos.total_adjusted_sealevel./1000,'b');
plot(roms.pos.time,roms.pos.zeta-mean(roms.pos.zeta,'omitnan'),'r');
plot(roms.pos.operational.time,roms.pos.operational.zeta-mean(roms.pos.operational.zeta,'omitnan'),'g');
% legend('raw SLX','adjusted SLX','ROMS-cwa','ROMS-cwa-OPERATIONAL','Location','best')
plot(gn.ts,gn.wl-.05,'k')
legend('adjusted SLX','ROMS-cwa','ROMS-cwa-OPERATIONAL','GER tide gauge','Location','best')
title('POS site Abrolhos')
ylabel('sealevel (m)')
set(gca,'fontsize',14)
set(gca,'xlim',[datenum(1992,1,1) datenum(2022,12,1)]);
datetick('x','yyyy','keeplimits')
set(gca,'ylim',[-1 1.2])
% datetick('x','yyyy/mm','keeplimits')
% 
% 
% figure; plot(slx.time,slx.total_adjusted_sealevel./1000);datetick
% hold on;
% plot(roms.time,roms.zeta-mean(roms.zeta),'r');
% raw.zeta=raw.z-mean(raw.z,'omitnan');
% 
% figure; plotraw.zeta
% figure; plot(raw.zeta)
% hol don
% plot(raw.time,raw.zeta,'g')
% legend('SLX','ROMS-cwa', 'raw SLX')
% datetick('x','yyyy-mmm','keeplimits')
% set(gca,'fontsize',14)
% grid on
% f=gcf
% f.color='w'
% f.Color='w'
% datetick('x','yyyy','keeplimits')
% figure; plot(slx.time,slx.total_adjusted_sealevel./1000);datetick
% hold on;
% plot(roms.time,roms.zeta-mean(roms.zeta),'r');
% set(gca,'fontsize',14)
% datetick('x','yyyy','keeplimits')
% legend('SLX','ROMS-cwa', 'raw SLX')
% f=gcf
% f.color='w'
% f.Color='w'
% datetick('x','yyyy/mm','keeplimits')


box on

output_filename=[whereput '/' site_id '_SLX_ROMS_sealevels']
savefig(output_filename)
 export_fig(output_filename,'-pdf','-painters');
%%  

load BOM_obs_archived-9.mat

north=archive(20)


figure; 
plot(north.mtime_UTC,north.press_msl)
datetick


%%  GERALDTON

% load SLX
ncref2='25665_Data.nc';
slx.ger.time=ncread(ncref2,'time')+datenum('1958-01-01');
slx.ger.total_adjusted_sealevel=double(ncread(ncref2,'sealevel'))./1000;

slx.ger.lat=ncread(ncref2,'latitude');
slx.ger.lon=ncread(ncref2,'longitude');
slx.ger.site_id=ncref2;
slx.ger.timezone=8;



% load roms reanalysis data
ncref='/Users/00068592/GitHub/cwa-roms/GER_2001-2020_qck_zeta.nc'
roms.ger.time=ncread(ncref,'ocean_time')./3600/24 +datenum(2000,1,1);
roms.ger.zeta=squeeze(ncread(ncref,'zeta'));

ncref='/Users/00068592/GitHub/cwa-roms/GER_2020-202205_OPERATIONAL_qck_zeta.nc'
roms.ger.operational.time=ncread(ncref,'ocean_time')./3600/24 +datenum(2000,1,1);
roms.ger.operational.zeta=squeeze(ncread(ncref,'zeta'));


%
% plot them
f=figure; 
f.Position=[69 439 2143 738];
f.Color='w';

% plot(raw.time,raw.zeta-mean(raw.zeta,'omitnan'),'g');
hold on;
plot(slx.ger.time,slx.ger.total_adjusted_sealevel,'b');
plot(roms.ger.time,roms.ger.zeta-mean(roms.ger.zeta,'omitnan'),'r');
plot(roms.ger.operational.time,roms.ger.operational.zeta-mean(roms.ger.operational.zeta,'omitnan'),'g');
plot(gn.ts,gn.wl-.05,'k')
% plot(slx.time,slx.total_adjusted_sealevel./1000,'r');
% legend('raw SLX','adjusted SLX','ROMS-cwa','ROMS-cwa-OPERATIONAL','Location','best')
legend('SLX','ROMS-cwa','ROMS-cwa-OPERATIONAL','tide gauge','pos slx','Location','best')
set(gca,'ylim',[-1 1.2])
title('GER')
ylabel('sealevel (m)')
set(gca,'fontsize',14)
set(gca,'xlim',[datenum(1992,1,1) datenum(2022,12,1)]);
datetick('x','yyyy','keeplimits')
box on
grid on

site_id=[slx.ger.site_id(1:5) '_GER'];
output_filename=[whereput '/' site_id '_SLX_ROMS_sealevels']
savefig(output_filename)
 export_fig(output_filename,'-pdf','-painters');
%% save sealevels to csv files

whereput=[pwd filesep 'SealevelX_ROMS_sealevels_v2']
mkdir(whereput);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             slx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[years months days hours minutes seconds]=datevec(slx.ger.time);
site_id=[slx.ger.site_id(1:5) '_GER'];
%  save to CSV comma delimted file
disp('Saving data to csv file...')
output_filename=[whereput '/' site_id '_sealevelX_sealevels']


    fid = fopen([output_filename '.csv'],'w');
%     
%     
    fprintf(fid,'Ivica SLX schism water levels for Jo Buckee study sites \n');
    fprintf(fid,'Created by Yasha Hetzel:               %s (UTC+8)\n',datestr(now));
    fprintf(fid,'With:                                  %s\n','plot_Jo_SLX_ROMS_water_levels.m');
    fprintf(fid,'Name:                                  %s\n','Post Office Island, Abrolhos WA');
    fprintf(fid,'ID:                                    %s\n',site_id);
    fprintf(fid,'Nearest grid point Longitude           %0.4f\n',slx.ger.lon);
    fprintf(fid,'Nearest grid point Latitude            %0.4f\n',slx.ger.lat);
    fprintf(fid,'Times are in UTC/GMT, To convert to local time add: %s hours\n',num2str(slx.ger.timezone));
    fprintf(fid,'\n');
    fprintf(fid,'Year, Month, Day, Hour, Total sea level (m)\n');
    fprintf(fid, '%d,%d,%d,%d,%0.4f\n',   [years, months, days, hours, slx.ger.total_adjusted_sealevel]')        
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       roms reanalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[years months days hours minutes seconds]=datevec(roms.ger.time);
site_id=[slx.ger.site_id(1:5) '_GER'];
%  save to CSV comma delimted file
disp('Saving data to csv file...')
output_filename=[whereput '/' site_id '_ROMS_REANALYSIS_sealevels']


    fid = fopen([output_filename '.csv'],'w');
%     
%     
    fprintf(fid,'Ivica Cwa-Roms REANALYSIS water levels for Jo Buckee study sites \n');
    fprintf(fid,'Created by Yasha Hetzel:               %s (UTC+8)\n',datestr(now));
    fprintf(fid,'With:                                  %s\n','plot_Jo_SLX_ROMS_water_levels.m');
    fprintf(fid,'Name:                                  %s\n','Geraldton WA');
    fprintf(fid,'ID:                                    %s\n',site_id);
    fprintf(fid,'Nearest grid point Longitude           %0.4f\n',slx.ger.lon);
    fprintf(fid,'Nearest grid point Latitude            %0.4f\n',slx.ger.lat);
    fprintf(fid,'Times are in UTC/GMT, To convert to local time add: %s hours\n',num2str(slx.ger.timezone));
    fprintf(fid,'\n');
    fprintf(fid,'Year, Month, Day, Hour, Total sea level (m)\n');
    fprintf(fid, '%d,%d,%d,%d,%0.4f\n',   [years, months, days, hours, roms.ger.zeta-mean(roms.ger.zeta,'omitnan')]')        
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       roms operational
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[years months days hours minutes seconds]=datevec(roms.ger.operational.time);
site_id=[slx.ger.site_id(1:5) '_GER'];
%  save to CSV comma delimted file
disp('Saving data to csv file...')
output_filename=[whereput '/' site_id '_ROMS_OPERATIONAL_sealevels']


    fid = fopen([output_filename '.csv'],'w');
%     
%     
    fprintf(fid,'Ivica Cwa-Roms OPERATIONAL water levels for Jo Buckee study sites \n');
    fprintf(fid,'Created by Yasha Hetzel:               %s (UTC+8)\n',datestr(now));
    fprintf(fid,'With:                                  %s\n','plot_Jo_SLX_ROMS_water_levels.m');
    fprintf(fid,'Name:                                  %s\n','Geraldton WA');
    fprintf(fid,'ID:                                    %s\n',site_id);
    fprintf(fid,'Nearest grid point Longitude           %0.4f\n',slx.ger.lon);
    fprintf(fid,'Nearest grid point Latitude            %0.4f\n',slx.ger.lat);
    fprintf(fid,'Times are in UTC/GMT, To convert to local time add: %s hours\n',num2str(slx.ger.timezone));
    fprintf(fid,'\n');
    fprintf(fid,'Year, Month, Day, Hour, Total sea level (m)\n');
    fprintf(fid, '%d,%d,%d,%d,%0.4f\n',   [years, months, days, hours, roms.ger.operational.zeta-mean(roms.ger.operational.zeta,'omitnan')]')        
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       tide gauge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[years months days hours minutes seconds]=datevec(gn.ts);
site_id=gn.name;
%  save to CSV comma delimted file
disp('Saving data to csv file...')
output_filename=[whereput '/' site_id '_TideGauge_sealevels']


    fid = fopen([output_filename '.csv'],'w');
%     
%     
    fprintf(fid,'Ivica Cwa-Roms OPERATIONAL water levels for Jo Buckee study sites \n');
    fprintf(fid,'Created by Yasha Hetzel:               %s (UTC+8)\n',datestr(now));
    fprintf(fid,'With:                                  %s\n','plot_Jo_SLX_ROMS_water_levels.m');
    fprintf(fid,'Name:                                  %s\n','Tide gauge observed Geraldton WA');
    fprintf(fid,'ID:                                    %s\n',site_id);
    fprintf(fid,'Nearest grid point Longitude           %0.4f\n',gn.lon);
    fprintf(fid,'Nearest grid point Latitude            %0.4f\n',gn.lat);
    fprintf(fid,'Times are in UTC/GMT, To convert to local time add: %s hours\n',num2str(slx.ger.timezone));
    fprintf(fid,'\n');
    fprintf(fid,'Year, Month, Day, Hour, Total sea level (m)\n');
    fprintf(fid, '%d,%d,%d,%d,%0.4f\n',   [years, months, days, hours, gn.wl]')        
fclose(fid);


%% 

output_filename=[whereput '/' 'Combined_POS_GER_sealevel_data'];
save(output_filename)




