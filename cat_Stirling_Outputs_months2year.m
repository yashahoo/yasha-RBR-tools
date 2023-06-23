% concatinate hourly outputs (monthly) from wave analysis or stirling
% created by process_Stirling_waves_sealevel_v9.m into year data.
clear all; 
close all;
clc;
addpath(genpath(pwd))

year = 2022 %2022;
location='Success';
datum_correction= [NaN] ;
README_time_adjustment = 'Time was corrected to local AWST by adding 8 hours'%'Time was corrected to local AWST  by adding 19 hours for may, June and 12 hours for other months'
Offset_from_UTC = '8';
instrument='FPA_pressure_sensor'
AnalysisMethod = 'zerocross' % 'zerocross' % 'spectral'
save_matfile       = 'Y'            % save the processed data to matfile
save_matfig        = 'N';           % save the .fig files?
print_fig          = 'Y';           % print to image file?
use_export_fig     = 'N'            % if print_fig ='Y' -->  'N' uses native matlab print function (default) OR 'Y' uses export_fig to make nice plots
img_type           = {'pdf','png'}; % if print_fig ='Y' -->  {'pdf','png','eps'} % options. can be multiple. pdf will not save if native matlab unless edit function
create_gapfree_time= 'N'            % to interpolate to perfect gap free time vector [ on raw data]- will cause problems with zero cross method if there are big gaps even if you fill them in!! %  NOT SURE HOW TO DEAL WITH THIS
max_timegap_days   = 5/60/24        % maximum time gap (in das) to interpolate data across when create_gapfree_time= 'Y'; 10/3600/24 = 10 seconds %  NOT SURE HOW TO DEAL WITH THIS

% whereput_combined_outputs=['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2020_Stirling_CLEANED_v4/' location '_' num2str(year) '_processed_combined_wave_data'];
whereput_combined_outputs=['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_' num2str(year) '_' location '_CLEANED_v1/' location '_' num2str(year) '_processed_combined_wave_data'];
mkdir(whereput_combined_outputs)


clear d dir
co=0;
for mo=1:12
    co=co+1;
%     zero-crossing
%     d(co)=dir(['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2020_Stirling_CLEANED_v4/2020' sprintf('%02.f',mo) '_Cleaned_pressure_data_Stirling_processed_' AnalysisMethod '/FPA_pressure_sensor_Stirling_PROCESSED_DATA_2020' sprintf('%02.f',mo) '*']);
% %     d(co)=dir(['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_' num2str(year) '_' location '_CLEANED_v1/' num2str(year) sprintf('%02.f',mo) '_' location '_Cleaned_pressure_data_' location '_processed_' AnalysisMethod '/FPA_pressure_sensor_' location '_PROCESSED_DATA_' num2str(year) sprintf('%02.f',mo) '*']);
% 
%     d(co)=dir(['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2022_Stirling_CLEANED_v2/2022' sprintf('%02.f',mo) '_Cleaned_pressure_data_Stirling_processed_' AnalysisMethod '/FPA_pressure_sensor_Stirling_PROCESSED_DATA_2022' sprintf('%02.f',mo) '*']);

    d(co)=dir(['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_' num2str(year-1) '-' num2str(year+1) '_' location '_CLEANED_v1/' num2str(year) sprintf('%02.f',mo) '_Cleaned_pressure_data_' location '_processed_' AnalysisMethod '/FPA_pressure_sensor_' location '_PROCESSED_DATA_' num2str(year) sprintf('%02.f',mo) '*']);
end
% dir();

m=matfile([d(1).folder filesep d(1).name]); % get first data to start concatinating
d(1).name
T_sl=m.T_sl; %Waves_output_Table
Waves_output_Table=m.Waves_output_Table;

for i=2:length(d)
    m=[];
    d(i).name
    m=matfile([d(i).folder filesep d(i).name]);
    tmp_sl=m.T_sl;
    tmp_Waves_output_Table=m.Waves_output_Table;
    
    T_sl=vertcat(T_sl,tmp_sl);
    
    Waves_output_Table=vertcat(Waves_output_Table,tmp_Waves_output_Table);

end






%% create QC flag for noisy IG events

remove_outliers='Y'     % replace outliers with nans? This uses moving window to identify outliers.... makes it much slower.
outlier_stds=4          % use this number of standard deviations from mean to identify bad data [10 is good...needs to be big to avoid throwing out real waves]
out_window = [3]  % moving window definition [behind infront] [600 600] is centred 10-minute window with 2 Hz data. If only single number, window is balanced
replace_nans='N'        % replace nans  by interpolation
% maxgap_secs=3700        % maximum gap in seconds to interpolate across [not used]


[IG_badi] = get_outliers_yh(Waves_output_Table.Hig,outlier_stds,out_window);
% bad_IG_flag=IG_badi;

 %% add the QC flag for bad IG waves to the table
 Waves_output_Table.Possible_bad_IG_flag=IG_badi;
 
% Waves_output_Table.bad_IG_flag=[];
% Waves_output_Table.Possible_bad_IG_flag=[];

%% get monthly x-ticks for long datasets

% wtt = table2timetable(Waves_output_Table,'RowTimes',datetime(Waves_output_Table.timestr));
% startDay = dateshift(wtt.Time(1), 'start', 'month');
% endDay =  dateshift(wtt.Time(end), 'end', 'month');


% or
wdt=datetime(Waves_output_Table.timestr);
startDay = dateshift(wdt(1), 'start', 'month');
endDay =  dateshift(wdt(end), 'end', 'month');
endYear =  dateshift(wdt(1), 'end', 'year');
Xtimes = [startDay:calmonths:endDay];
[y1,m1,d1]=ymd(startDay);[y2,m2,d2]=ymd(endYear); 
Xlims= [startDay endYear];




%% plot IG waves seperately

IGbi=find(Waves_output_Table.Possible_bad_IG_flag>0);
clean_IG=Waves_output_Table.Hig;
clean_IG(IGbi)=NaN;




  f=figure;
 f.Position=[119 661 2082 649];
 f.Color='w'
 plot(datenum(Waves_output_Table.timestr),Waves_output_Table.Hig);
 hold on;
  plot(datenum(Waves_output_Table.timestr(IGbi,:)),Waves_output_Table.Hig(IGbi),'r*')
%  title([num2str(m.subsample_mins) ' minute averaged water level'])
ylabel('Infragravity (30-200 secs) height (m)')
set(gca,'fontsize',14)
xlabel(['(' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy') ')']);
legend('IG band waves','Possible_bad_IG_flag','interpreter','none')
 ax=gca;
 ax.XTick=datenum(Xtimes);
 datetick('x','mmm','keeplimits','keepticks');
fname=[whereput_combined_outputs '/' location ' _' AnalysisMethod '_IGwaves_flags_hourly_' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy')];
% print('-dpng','-r200', fname);
% print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)

 
   f=figure;
 f.Position=[119 661 2082 649];
 f.Color='w'
 plot(datenum(Waves_output_Table.timestr),clean_IG);
 hold on;
%   plot(datenum(Waves_output_Table.timestr(IGbi,:)),Waves_output_Table.Hig(IGbi),'r*')
%  title([num2str(m.subsample_mins) ' minute averaged water level'])
ylabel('Cleaned Infragravity (30-200 secs) height (m)')
 set(gca,'fontsize',14)
 xlabel(['(' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy') ')']);
  ax=gca;
 ax.XTick=datenum(Xtimes);
 datetick('x','mmm','keeplimits','keepticks');
 

 fname=[whereput_combined_outputs '/' location ' _' AnalysisMethod '_IGwaves_CLEANED_hourly_' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy')];
% print('-dpng','-r200', fname);
% print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
 
 



%% adjust vertical datum for feb, apr,august

%  For 1 min water levels
% 
% tstr1=T_sl.timestr;
% ts1=datenum(tstr1);
% p1=T_sl.sub_press;
% febi=find(ts1>datenum(2022,1,31,23,59,59) & ts1<datenum(2022,3,1));
% 
% t2=mean(p1(febi(1:7200)),'omitnan')
% t1=mean(p1(febi(1)-7200:febi(1)-1),'omitnan')
% dz1=t2-t1
% 
% p1(febi)=p1(febi)-dz1;
% 
% figure; plot(ts1,p1);datetick
% 
% apri=find(ts1>datenum(2022,3,31,23,59,59) & ts1<datenum(2022,5,1));
% t2=mean(p1(apri(1:7200)),'omitnan')
% t1=mean(p1(apri(1)-7200:apri(1)-1),'omitnan')
% dz1=t2-t1
% 
% p1(apri)=p1(apri)-dz1; %+0.02 %
% figure; plot(ts1,p1);datetick
% 
% 
% augi=find(ts1>datenum(2022,7,31,23,59,59) & ts1<datenum(2022,9,1));
% t2=mean(p1(augi(1:7200)),'omitnan')
% t1=mean(p1(augi(1)-7200:augi(1)-1),'omitnan')
% dz1=t2-t1
% 
% p1(augi)=p1(augi)-dz1; %+0.02 %
% figure; plot(ts1,p1);datetick
% 
% % T_sl.sub_press=p1;

%% adjust vertical datum for feb, apr,august for 1 hr water levels for 2022
% 
% tstr1=Waves_output_Table.timestr;
% ts1=datenum(tstr1);
% p1=Waves_output_Table.Sealevel;
% febi=find(ts1>datenum(2022,1,31,23,59,59) & ts1<datenum(2022,3,1));
% 
% t2=mean(p1(febi(1:120)),'omitnan')
% t1=mean(p1(febi(1)-120:febi(1)-1),'omitnan')
% dz1=t2-t1
% 
% p1(febi)=p1(febi)-dz1;
% 
% figure; plot(ts1,p1);datetick
% 
% apri=find(ts1>datenum(2022,3,31,23,59,59) & ts1<datenum(2022,5,1));
% t2=mean(p1(apri(1:120)),'omitnan')
% t1=mean(p1(apri(1)-120:apri(1)-1),'omitnan')
% dz1=t2-t1
% 
% p1(apri)=p1(apri)-dz1; %+0.02 %
% figure; plot(ts1,p1);datetick
% 
% 
% augi=find(ts1>datenum(2022,7,31,23,59,59) & ts1<datenum(2022,9,1));
% t2=mean(p1(augi(1:120)),'omitnan')
% t1=mean(p1(augi(1)-120:augi(1)-1),'omitnan')
% dz1=t2-t1
% 
% p1(augi)=p1(augi)-dz1; %+0.02 %
% figure; plot(ts1,p1);datetick
% 
% 
% % Waves_output_Table.Sealevel=p1;






 %% save matfile and write the tables to text files
 
readme=['Hourly wave outputs saved as matlab table: Waves_output_Table; subsampled sealevel (pressure,tide,residual)saved in matlab table: T_sl. Timezone = UTC+' Offset_from_UTC]
outf=[whereput_combined_outputs '/' location '_waves_hourly_' AnalysisMethod '_' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy')];
save(outf,'Waves_output_Table','T_sl','readme','README_time_adjustment')
                  
outf_tname=[whereput_combined_outputs '/' location '_waves_hourly_' AnalysisMethod '_' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy') '.csv'];
writetable(Waves_output_Table,outf_tname,'Delimiter',',','QuoteStrings',true)
disp(['Hourly wave and sea level data saved to csv text file: ' outf_tname])

% dtminstr=num2str(mean(diff(datenum(T_sl.timestr(1:5,:))))*24*60);
dtminstr=num2str(round(mean(diff(datenum(T_sl.timestr(1:5,:))))*24*60));
outf_tname=[whereput_combined_outputs '/' location '_' dtminstr '_min_' 'SEALEVEL_'  datestr(mean(datenum(T_sl.timestr(1,:)),'omitnan'),'yyyy') '.csv'];
writetable(T_sl,outf_tname,'Delimiter',',','QuoteStrings',true)
disp([dtminstr ' minute wave and sea level data saved to csv text file: ' outf_tname])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   interpolate 1 min wate rlevels to perfect time vector and then save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % %           do for sub-hourly data
% %             
% %             first_good_ti=find(~isnan(ts),1,'first');
% %             last_good_ti=find(~isnan(ts),1,'last');
% %             t_interp=[ts(first_good_ti):datenum(0,0,0,0,0,1/sample_rate_hz):ts(last_good_ti)];
% %             p_interp=interp1gap(ts,p,t_interp,max_timegap_days,'linear','extrap',NaN); % 'nearest'
% %             
% %             ts=t_interp';
% %             p=p_interp';
% %             
% %             disp('INTERPOLATING TO GAP FREE TIME -- BEWARE!!!!')


% wtt = table2timetable(Waves_output_Table,'RowTimes',datetime(Waves_output_Table.timestr));
% startDay = dateshift(wtt.Time(1), 'start', 'day');
% endDay =  dateshift(wtt.Time(end), 'end', 'day');
% newTimes = [startDay:hours(1):endDay];
% 
% 
% wtt_backup=wtt;
% newtime = wtt.Time+minutes(13);
% wtt.Time=newtime;
% 
% % wtt.Properties.VariableContinuity = {'unset','continuous','continuous','continuous','continuous','continuous','continuous','continuous','continuous','continuous','continuous','continuous','step'};
% wtt.Properties.VariableContinuity = {'unset','unset','unset','unset','unset','unset','unset','unset','unset','unset','unset','unset','step'};
% wtt.Properties
% 
% ttr=retime(wtt,newTimes,'nearest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   interpolate to perfect HOURLY time vector and then save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% % this is tricker than i thought bc of the gap size... i want to only interpolateto a tiem vector if teh gap not too big
% 
% wtt = table2timetable(Waves_output_Table,'RowTimes',datetime(Waves_output_Table.timestr));
% startDay = dateshift(wtt.Time(1), 'start', 'day');
% endDay =  dateshift(wtt.Time(end), 'end', 'day');
% newTimes = [startDay:hours(1):endDay];
% 
% 
% wtt_backup=wtt;
% newtime = wtt.Time+minutes(13);
% wtt.Time=newtime;
% 
% % wtt.Properties.VariableContinuity = {'unset','continuous','continuous','continuous','continuous','continuous','continuous','continuous','continuous','continuous','continuous','continuous','step'};
% wtt.Properties.VariableContinuity = {'unset','unset','unset','unset','unset','unset','unset','unset','unset','unset','unset','unset','step'};
% wtt.Properties
% 
% ttr=retime(wtt,newTimes,'nearest');
% 
% % 
% %  fill except where gaps
% % % Find gaps and create group vector
% idx = diff(wtt.Time) > hours(1);
% idx = [true; idx];
% group = cumsum(idx);
% % % Apply retime function for each group and concatenate the result
% TT2 = [];
% for kk = 1:max(group)
%   idx = group == kk;
%   TT2 = [TT2; retime(wtt(idx,:),newTimes,'default')];
% end


% 
% % subset  species for processing
% ispecies=strcmpi(species,w.Species);
% 
% %  create gap-free time
% sw=w(ispecies,:); % species specific
% sw.catchtime=datetime(sw.Year,sw.Month_Caught,15);
% sw.spawntime=datetime(sw.Year,sw.Month_Spawning,15);
% 
% % get only the spawning times
% coli=(strcmpi(sw.Properties.VariableNames,'NSB_Spawning_weight') | strcmpi(sw.Properties.VariableNames,'DS_Spawning_weight') | strcmpi(sw.Properties.VariableNames,'Species'))
% spawn = table2timetable(sw(:,coli),'RowTimes',sw.spawntime);
% spawntime_complete = [datetime(spawn.Time.Year(1),1,15):calmonths(1):datetime(spawn.Time.Year(end),12,15)]';
% spawni=retime(spawn,spawntime_complete);
% 

% 
% % for rottnest need to convert the timetable back to table as in r2015a not
% % supported
% King_spawni=timetable2table(spawni);
% Tiger_spawni=timetable2table(spawni);
% save('Spawning_weight_tables','King_spawni','Tiger_spawni')






%     switch create_gapfree_time 
%         case 'Y'
%             
%             switch AnalysisMethod
%             
%                 case 'zerocross'
%             
%                  ts=datenum(   
%                     
%                  
%                 case 'spectral'
%                     
%                     
%             end   
%                  
%             
% %           do for sub-hourly data
%             
%             first_good_ti=find(~isnan(ts),1,'first');
%             last_good_ti=find(~isnan(ts),1,'last');
%             t_interp=[ts(first_good_ti):datenum(0,0,0,0,0,1/sample_rate_hz):ts(last_good_ti)];
%             p_interp=interp1gap(ts,p,t_interp,max_timegap_days,'linear','extrap',NaN); % 'nearest'
%             
%             ts=t_interp';
%             p=p_interp';
%             
%             disp('INTERPOLATING TO GAP FREE TIME -- BEWARE!!!!')
%             
%     end
    


%% plot water level 
ts=datenum(Waves_output_Table.timestr);
Htime=ts;
burst_times=Htime;

% x1=datenum(str2num(datestr(min(ts),'yyyy')),str2num(datestr(min(ts),'mm')),1);
% x2=datenum(str2num(datestr(min(ts),'yyyy')),str2num(datestr(max(ts),'mm')),31,23,59,59);
% xlims=[x1 x2]; % limits
% xt=datenum(str2num(datestr(min(ts),'yyyy')),1:12,1); % tick marks




%  f=figure;
%  f.Position=[119 661 2082 649];
%  f.Color='w'
%  plot(datenum(Waves_output_Table.timestr),Waves_output_Table.Sealevel);datetick
%  title([num2str(diff(datenum(Waves_output_Table.timestr(1:2,:)))*24*60) ' minute averaged water level'])
%  set(gca,'fontsize',14)
%  xlabel(['(' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy') ')']);
 
 
 time_subsampled=datenum(T_sl.timestr);
 press_subsampled=T_sl.sub_press;
 tide_subsampled=T_sl.sub_tide;
 residual_subsampled=T_sl.sub_residual;
 dud=press_subsampled-nanmean(press_subsampled);
 
 f=figure;
 f.Color='w';
 f.Position=[119 661 2082 649];
 hold on
 plot(time_subsampled,press_subsampled-nanmean(press_subsampled),'k')
 plot(time_subsampled,tide_subsampled,'b')
 plot(time_subsampled,residual_subsampled,'r')
 plot([time_subsampled(1) time_subsampled(end)],[0 0],'color',[.5 .5 .5],'linewidth',0.5)
 legend('total','tide','residual','location','best')
 ylabel('Water level (m)')
 %         xlabel(['Date (' datestr(time_subsampled(1),'YYYY') ')'])
 yt=[floor(min(dud)):.25:ceil(max(dud))];
 set(gca,'ytick',yt)
ax=gca;
 ax.XTick=datenum(Xtimes);
%  ax.XLim=[datenum(startDay) datenum(endDay)];
  ax.XLim=[datenum(startDay) datenum(endYear)];
 datetick('x','mmm','keeplimits','keepticks');
%  datetick('x','dd-mmm','keeplimits','keepticks');
 ax.XTickLabelRotation=90;
 box on;
 grid on;
 set(gca,'fontsize',14);
 xlabel(datestr(mean(time_subsampled),'YYYY'))
 
 fname=[whereput_combined_outputs '/' location '_SEALEVEL_ONLY_' dtminstr '_min_' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy')];
 % print('-dpng','-r200', fname);
 % print('-dpdf', '-painters', fname);
 %         export_fig(fname,'-pdf','-png','-transparent')
 savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
 
 %% plot wave height alone
 
 ts=datenum(Waves_output_Table.timestr);
Htime=ts;
burst_times=Htime;

x1=datenum(str2num(datestr(min(ts),'yyyy')),str2num(datestr(min(ts),'mm')),1);
x2=datenum(str2num(datestr(min(ts),'yyyy')),str2num(datestr(max(ts),'mm')),31,23,59,59);
xlims=[x1 x2]; % limits
xt=datenum(str2num(datestr(min(ts),'yyyy')),1:12,1); % tick marks


 switch AnalysisMethod
     case 'zerocross'
 f=figure;
 f.Position=[119 661 2082 649];
 f.Color='w'
 plot(datenum(Waves_output_Table.timestr)+0/24,Waves_output_Table.Hs);
 hold on
         legend('Hs')
         ylabel('Wave height (m)')
         set(gca,'xlim',xlims,'xtick',xt)
         set(gca,'ylim',[0 ceil(max(Waves_output_Table.Hs))])
         ax=gca;
 ax.XTick=datenum(Xtimes);
%  ax.XLim=[datenum(startDay) datenum(endDay)];
  ax.XLim=[datenum(startDay) datenum(endYear)];
 datetick('x','mmm','keeplimits','keepticks');
%  datetick('x','dd-mmm','keeplimits','keepticks');
 ax.XTickLabelRotation=90;
 box on;
 grid on;
 set(gca,'fontsize',14);
 xlabel(datestr(mean(time_subsampled),'YYYY'))
%  title([num2str(m.subsample_mins) ' minute averaged water level'])
% ylabel('Hs (m)')
 set(gca,'fontsize',14)
 xlabel(['(' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy') ')']);
                 fname=[whereput_combined_outputs '/' location ' _' AnalysisMethod '_WAVES_ONLY_hourly_' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy')];
                % print('-dpng','-r200', fname);
                % print('-dpdf', '-painters', fname);
                %         export_fig(fname,'-pdf','-png','-transparent')
                savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)      
 
 
     case 'spectral'
         
         f=figure;
         f.Position=[119 661 2082 649];
         f.Color='w'
         plot(burst_times,Waves_output_Table.Hm0swell)
         hold on
         plot(ts,Waves_output_Table.Hm0sea)
         plot(ts,Waves_output_Table.Hm0,'k')
         legend('Hm0swell','Hm0sea','Hm0')
         ylabel('Wave height (m)')
%          set(gca,'xlim',xlims,'xtick',xt)
         set(gca,'ylim',[0 ceil(max(Waves_output_Table.Hm0))])
         ax=gca;
         ax.XTick=datenum(Xtimes);
         %  ax.XLim=[datenum(startDay) datenum(endDay)];
         ax.XLim=[datenum(startDay) datenum(endYear)];
         datetick('x','mmm','keeplimits','keepticks');
         %  datetick('x','dd-mmm','keeplimits','keepticks');
         ax.XTickLabelRotation=90;
         box on;
         grid on;
         set(gca,'fontsize',14);
         xlabel(datestr(mean(time_subsampled),'YYYY'))
         ax.XTickLabelRotation=90;
         set(gca,'fontsize',14)
         %                             xlabel(datestr(mean(ts),'YYYY'))
         title([instrument ' ' location ],'interpreter','none')
         box on;
         grid on;
         xlabel(['(' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy') ')']);
%          ax2=subaxis(4,1,2)
%          %                             h1=plot(ts,Waves_output_Table.Tp)
%          hold on
%          h2=plot(ts,Waves_output_Table.Tpswell)
%          h3=plot(ts,Waves_output_Table.Tpsea)
%          set(gca,'xlim',xlims,'xtick',xt,'ytick',[0:2:24],'ylim',[0 22])
%          %                             legend('Tp','swell','sea');
%          legend('swell','sea');
%          ylabel('Period (s)')
%          %                             title([instrument ' ' location ],'interpreter','none')
%          datetick('x','dd-mmm','keeplimits','keepticks')
%          %         xlabel(datestr(mean(ts),'YYYY'))
%          set(gca,'xticklabel',[])
%          %                             if exist('Offset_from_UTC','var')
%          %                                 xlabel([datestr(mean(ts),'YYYY') ' (UTC+' Offset_from_UTC ')'])
%          %                             else
%          %                                 xlabel(datestr(mean(ts),'YYYY'))
%          %                             end
%          %         rotateXLabels( gca(), 90 )
%          ax2.XTickLabelRotation=90;
%          %                             h1.Color=[0 0 0];
%          h2.Color=[ 0    0.4470    0.7410];
%          h3.Color=[0.8500    0.3250    0.0980];
         %                             box on;
         %                             grid on;
         
         
                fname=[whereput_combined_outputs '/' location ' _' AnalysisMethod '_WAVES_ONLY_hourly_' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy')];
                % print('-dpng','-r200', fname);
                % print('-dpdf', '-painters', fname);
                %         export_fig(fname,'-pdf','-png','-transparent')
                savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)       
         
         
 end
 
 % % %  
% % %  %  plot IG
% % %  f=figure;
% % %  f.Position=[119 661 2082 649];
% % %  f.Color='w'
% % %  plot(datenum(Waves_output_Table.timestr),Waves_output_Table.Hig);
% % %  datetick
% % %  %  title([num2str(m.subsample_mins) ' minute averaged water level'])
% % %  ylabel('Infragravity (30-200 secs) height (m)')
% % %  set(gca,'fontsize',14)
% % %  xlabel(['(' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy') ')']);
% % %  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot all [ need to edit for wave outputs]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts=datenum(Waves_output_Table.timestr);
Htime=ts;
burst_times=Htime;
time_subsampled=datenum(T_sl.timestr);
press_subsampled=T_sl.sub_press;
tide_subsampled=T_sl.sub_tide;
residual_subsampled=T_sl.sub_residual;


% 
% x1=datenum(str2num(datestr(min(ts),'yyyy')),str2num(datestr(min(ts),'mm')),1);
% x2=datenum(str2num(datestr(min(ts),'yyyy')),str2num(datestr(max(ts),'mm')),31,23,59,59);
% xlims=[x1 x2]; % limits
% xt=datenum(str2num(datestr(min(ts),'yyyy')),1:12,1); % tick marks


calculate_waves = 'Y'
% % % AnalysisMethod = 'zerocross' %'spectral'
SeparateSeaSwell = 'yes' % only for spectral
datum = 'MSL' %'CD'
begin=min(ts);
finish=max(ts);


    switch calculate_waves
        case 'Y'
            % set  date ranges
%             xlim(1,:)=[begin finish]; % all
            % xlim(2,:)=[datenum(2017,9,4)  datenum(2017,9,10)];
            % xlim(3,:)=[datenum(2017,9,13)  datenum(2017,9,16)];
            % xlim(4,:)=[datenum(2017,10,1)  datenum(2017,10,25)];
            % xlim(5,:)=[datenum(2017,10,13)  datenum(2017,10,16)];
            
%             clear xlims
            
            for i=1 %:5
%                 xlims=xlim(i,:);
%                 xt=[xlims(1):xlims(2)];
                
                
                f=figure;
                f.Position=[200 100 600 800];
                f.Color ='w';
                
                switch AnalysisMethod
                    case 'spectral'
                        
                        switch SeparateSeaSwell
                            case 'yes'
                                
                                
                                ax1=subaxis(4,1,1)
                                plot(burst_times,Waves_output_Table.Hm0swell)
                                hold on
                                plot(ts,Waves_output_Table.Hm0sea)
                                plot(ts,Waves_output_Table.Hm0,'k')
                                legend('Hm0swell','Hm0sea','Hm0')
                                ylabel('Wave height (m)')
                                set(gca,'xlim',xlims,'xtick',xt)
                                set(gca,'ylim',[0 ceil(max(Waves_output_Table.Hm0))])
                                ax=gca;
 ax.XTick=datenum(Xtimes);
%  ax.XLim=[datenum(startDay) datenum(endDay)];
  ax.XLim=[datenum(startDay) datenum(endYear)];
 datetick('x','mmm','keeplimits','keepticks');
%  datetick('x','dd-mmm','keeplimits','keepticks');
 ax.XTickLabelRotation=90;
 box on;
 grid on;
%  set(gca,'fontsize',14);
 xlabel(datestr(mean(time_subsampled),'YYYY'))
                                set(gca,'xticklabel',[])
                                %         rotateXLabels( gca(), 90 )
                                ax.XTickLabelRotation=90;
                                %                             xlabel(datestr(mean(ts),'YYYY'))
                                title([instrument ' ' location ],'interpreter','none')
                                box on;
                                grid on;
                                
                                ax2=subaxis(4,1,2)
                                %                             h1=plot(ts,Waves_output_Table.Tp)
                                hold on
                                h2=plot(ts,Waves_output_Table.Tpswell)
                                h3=plot(ts,Waves_output_Table.Tpsea)
                                set(gca,'ytick',[0:2:24],'ylim',[0 22])
                                %                             legend('Tp','swell','sea');
                                legend('swell','sea');
                                ylabel('Period (s)')
                                %                             title([instrument ' ' location ],'interpreter','none')
%                                 datetick('x','dd-mmm','keeplimits','keepticks')
ax=gca;
 ax.XTick=datenum(Xtimes);
%  ax.XLim=[datenum(startDay) datenum(endDay)];
  ax.XLim=[datenum(startDay) datenum(endYear)];
 datetick('x','mmm','keeplimits','keepticks');
%  datetick('x','dd-mmm','keeplimits','keepticks');
 ax.XTickLabelRotation=90;
 box on;
 grid on;
%  set(gca,'fontsize',14);
 xlabel(datestr(mean(time_subsampled),'YYYY'))

                                %         xlabel(datestr(mean(ts),'YYYY'))
                                set(gca,'xticklabel',[])
                                %                             if exist('Offset_from_UTC','var')
                                %                                 xlabel([datestr(mean(ts),'YYYY') ' (UTC+' Offset_from_UTC ')'])
                                %                             else
                                %                                 xlabel(datestr(mean(ts),'YYYY'))
                                %                             end
                                %         rotateXLabels( gca(), 90 )
                                ax2.XTickLabelRotation=90;
                                %                             h1.Color=[0 0 0];
                                h2.Color=[ 0    0.4470    0.7410];
                                h3.Color=[0.8500    0.3250    0.0980];
                                %                             box on;
                                %                             grid on;
                                
                            case 'no'
                                
                                ax1=subaxis(4,1,1)
                                plot(ts,Waves_output_Table.Hm0,'k')
                                legend('Hm0')
                                ylabel('Wave height (m)')
                                set(gca,'xlim',xlims,'xtick',xt)
                                set(gca,'ylim',[0 ceil(max(Waves_output_Table.Hm0))])
%                                 datetick('x','dd-mmm','keeplimits','keepticks')
ax=gca;
 ax.XTick=datenum(Xtimes);
%  ax.XLim=[datenum(startDay) datenum(endDay)];
  ax.XLim=[datenum(startDay) datenum(endYear)];
 datetick('x','mmm','keeplimits','keepticks');
%  datetick('x','dd-mmm','keeplimits','keepticks');
 ax.XTickLabelRotation=90;
 box on;
 grid on;
%  set(gca,'fontsize',14);
 xlabel(datestr(mean(time_subsampled),'YYYY'))
                                set(gca,'xticklabel',[])
                                %         rotateXLabels( gca(), 90 )
                                ax1.XTickLabelRotation=90;
                                %                             xlabel(datestr(mean(ts),'YYYY'))
                                title([instrument ' ' location ],'interpreter','none')
                                box on;
                                grid on;
                                
                                ax2=subaxis(4,1,2)
                                h1=plot(ts,Waves_output_Table.Tp)
                                set(gca,'xlim',xlims,'xtick',xt,'ytick',[0:2:24],'ylim',[0 22])
                                legend('Tp','swell','sea');
                                ylabel('Period (s)')
                                %                             title([instrument ' ' location ],'interpreter','none')
                                datetick('x','dd-mmm','keeplimits','keepticks')
                                set(gca,'xticklabel',[])
                                %         xlabel(datestr(mean(ts),'YYYY'))
                                %                             if exist('Offset_from_UTC','var')
                                %                                 xlabel([datestr(mean(ts),'YYYY') ' (UTC+' Offset_from_UTC ')'])
                                %                             else
                                %                                 xlabel(datestr(mean(ts),'YYYY'))
                                %                             end
                                %         rotateXLabels( gca(), 90 )
                                ax2.XTickLabelRotation=90;
                                h1.Color=[0 0 0];
                                %                 h2.Color=[ 0    0.4470    0.7410];
                                %                 h3.Color=[0.8500    0.3250    0.0980];
                                %                             box on;
                                %                             grid on;
                                %
                                
                        end
                        
                        
                    case 'zerocross'
                        %              Hs: [100×1 double]
                        %              Hz: [100×1 double]
                        %              Tz: [100×1 double]
                        %              Ts: [100×1 double]
                        % Hmax
                        
                        ax1=subaxis(4,1,1)
                        plot(ts,Waves_output_Table.Hs)
                        hold on
                        plot(ts,Waves_output_Table.Hz)
%                         plot(ts,Waves_output_Table.Hmax)
                        legend('Hs','Hz','Hmax')
                        ylabel('Wave height (m)')
                        set(gca,'xlim',xlims,'xtick',xt,'ytick',[0:0.5:5])
                        set(gca,'ylim',[0 4]) 
                        ax1.YLim=[0 ceil(max(Waves_output_Table.Hs))]
%                         datetick('x','dd-mmm','keeplimits','keepticks')
ax=gca;
 ax.XTick=datenum(Xtimes);
%  ax.XLim=[datenum(startDay) datenum(endDay)];
  ax.XLim=[datenum(startDay) datenum(endYear)];
 datetick('x','mmm','keeplimits','keepticks');
%  datetick('x','dd-mmm','keeplimits','keepticks');
 ax.XTickLabelRotation=90;
 box on;
 grid on;
%  set(gca,'fontsize',14);
 xlabel(datestr(mean(time_subsampled),'YYYY'))
                        %         rotateXLabels( gca(), 90 )
                        ax1.XTickLabelRotation=90;
                        set(gca,'xticklabel',[])
                        %                     xlabel(datestr(mean(ts),'YYYY'))
                        title([instrument ' ' location ],'interpreter','none')
                        box on;
                        grid on;
                        
                        ax2=subaxis(4,1,2)
                        plot(ts,Waves_output_Table.Ts)
                        hold on
                        plot(ts,Waves_output_Table.Tz)
                        set(gca,'xlim',xlims,'xtick',xt,'ytick',[0:2:24],'ylim',[0 22])
                        legend('Ts','Tz')
                        ylabel('Period (s)')
                        %                     title([instrument ' ' location ],'interpreter','none')
ax2.YGrid='on'
ax2.XGrid='on'
ax2.XTick=datenum(Xtimes);
ax2.XLim=[datenum(startDay) datenum(endYear)];
datetick('x','mmm','keeplimits','keepticks');
                        set(gca,'xticklabel',[])
                        %                     if exist('Offset_from_UTC','var')
                        %                         xlabel([datestr(mean(ts),'YYYY') ' (UTC+' Offset_from_UTC ')'])
                        %                     else
                        %                         xlabel(datestr(mean(ts),'YYYY'))
                        %                     end
                        %         rotateXLabels( gca(), 90 )
                        ax2.XTickLabelRotation=90;
                        box on;
                        %                                         grid on;
                end
                % ---------------------------------------------------------------------- %
                % plot water levels
                % ---------------------------------------------------------------------- %
                
                switch datum
                    case 'CD'
                        % plot chart datum
                        ax3=subaxis(4,1,3);
                        hold on
                        plot(time_subsampled,press_subsampled,'k')
                        plot(time_subsampled,tide_subsampled+mean(press_subsampled),'b')
                        plot(time_subsampled,residual_subsampled,'r')
                        plot([time_subsampled(1) time_subsampled(end)],[0 0],'color',[.5 .5 .5],'linewidth',0.5)
                        % plot(time_subsampled(residual_subsampled<0),residual_subsampled(residual_subsampled<0),'g')
                        legend('total','tide','residual','location','best')
                        ylabel('Water level (m CD)')
                        %         xlabel(['Date (' datestr(time_subsampled(1),'YYYY') ')'])
                        yt=[floor(min(residual_subsampled)):.25:ceil(max(press_subsampled))];
                        set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
%                         datetick('x','dd-mmm','keeplimits','keepticks')
ax=gca;
 ax.XTick=datenum(Xtimes);
%  ax.XLim=[datenum(startDay) datenum(endDay)];
  ax.XLim=[datenum(startDay) datenum(endYear)];
 datetick('x','mmm','keeplimits','keepticks');
%  datetick('x','dd-mmm','keeplimits','keepticks');
 ax.XTickLabelRotation=90;
 box on;
 grid on;
%  set(gca,'fontsize',14);
 xlabel(datestr(mean(time_subsampled),'YYYY'))
                        set(gca,'xticklabel',[])
                        %         rotateXLabels( gca(), 90 )
                        ax3.XTickLabelRotation=90;
                        box on;
                        grid on;
                        
                    case 'MSL'
                        dud=press_subsampled-nanmean(press_subsampled);
                        % plot msl
                        ax3=subaxis(4,1,3);
                        hold on
                        plot(time_subsampled,press_subsampled-nanmean(press_subsampled),'k')
                        plot(time_subsampled,tide_subsampled,'b')
                        plot(time_subsampled,residual_subsampled,'r')
                        plot([time_subsampled(1) time_subsampled(end)],[0 0],'color',[.5 .5 .5],'linewidth',0.5)
                        legend('total','tide','residual','location','best')
                        ylabel('Water level (m)')
                        %         xlabel(['Date (' datestr(time_subsampled(1),'YYYY') ')'])
                        yt=[floor(min(dud)):.25:ceil(max(dud))];
%                         set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
ax=gca;
 ax.XTick=datenum(Xtimes);
%  ax.XLim=[datenum(startDay) datenum(endDay)];
  ax.XLim=[datenum(startDay) datenum(endYear)];
 datetick('x','mmm','keeplimits','keepticks');
%  datetick('x','dd-mmm','keeplimits','keepticks');
 ax.XTickLabelRotation=90;
 box on;
 grid on;
%  set(gca,'fontsize',14);
 xlabel(datestr(mean(time_subsampled),'YYYY'))
%                         datetick('x','dd-mmm','keeplimits','keepticks')
                        set(gca,'xticklabel',[])
                        %         rotateXLabels( gca(), 90 )
                        ax3.XTickLabelRotation=90;
                        box on;
                        grid on;
                        
                        %                 h2.Color=[ 0    0.4470    0.7410];
                        %                 h3.Color=[0.8500    0.3250    0.0980];
                        
                end
                
                
                
                
                % ---------------------------------------------------------------------- %
                
                
                % ---------------------------------------------------------------------- %
                %   plot IG waves
                % ---------------------------------------------------------------------- %
                ax4=subaxis(4,1,4)
                hold on;
%                 plot(Htime,Waves_output_Table.Hig)
                plot(Htime,clean_IG)
                plot([Htime(1) Htime(end)],[0 0],'color',[.5 .5 .5],'linewidth',0.5)
                legend('IG waves - outliers removed','location','best')
                set(gca,'xlim',xlims,'xtick',xt)
                set(gca,'ytick',[-1:.05:1],'ylim',[0 .3])
%                 datetick('x','dd-mmm','keeplimits','keepticks')
ax=gca;
 ax.XTick=datenum(Xtimes);
%  ax.XLim=[datenum(startDay) datenum(endDay)];
  ax.XLim=[datenum(startDay) datenum(endYear)];
 datetick('x','mmm','keeplimits','keepticks');
%  datetick('x','dd-mmm','keeplimits','keepticks');
 ax.XTickLabelRotation=90;
 box on;
 grid on;
%  set(gca,'fontsize',14);
 xlabel(datestr(mean(time_subsampled),'YYYY'))
                %             rotateXLabels( gca(), 90 )
                ax3.XTickLabelRotation=90;
                ylabel('Height (m)')
                
                if exist('Offset_from_UTC','var')
                    xlabel([datestr(mean(Htime),'YYYY') ' (UTC+' Offset_from_UTC ')'])
                else
                    xlabel(datestr(mean(Htime),'YYYY'))
                end
                
                grid on
                
ax1.XLabel  =[];
ax2.XLabel  =[];
ax3.XLabel  =[];
%                 fname=[whereput '/' instrument '_' location '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
                fname=[whereput_combined_outputs '/' location ' _' AnalysisMethod '_waves_sealevel_hourly_' datestr(mean(datenum(Waves_output_Table.timestr(1,:)),'omitnan'),'yyyy')];
                % print('-dpng','-r200', fname);
                % print('-dpdf', '-painters', fname);
                %         export_fig(fname,'-pdf','-png','-transparent')
                savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
                
            end
            
        case 'N'
            disp('NO waves plotted')
    end
    
    

    