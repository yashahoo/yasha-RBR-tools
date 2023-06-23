% divide chari's annual 2hz data from stirling to monthly matfiles
% clean it too
% yasha hetzel 2023-04-20

clear all;
clc;
close all;


restoredefaultpath % this just fixes if there are any conflicting functions
% % should work if you run it from the yasha_RBR_tools_v9 directory
addpath(genpath(pwd))







%%


%%-------------------------------------------------------------------------%
%%                          Settings
%%-------------------------------------------------------------------------%

years2do = [2021] % [2017 2018] % can be vector of years (for loop)
months2do = [1:2] %[7 8]   % can be vector of months (e.g. [1:12] for loop)

location = 'Parmelia'
% input_directory='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2020_Stirling'
% output_directory=['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2020_Stirling_CLEANED_v4']
input_directory='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2020_Parmelia'
output_directory=['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2020_Parmelia_CLEANED_v1']

% move problematic files to here:
% error_folder='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2020_Stirling/error_folder'
error_folder='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2020_Parmelia/error_folder'

sRateHz=2; % sampling rate in Hz
UTC_offset=8; % time in local AWST UTC+8hrs
fix_offset = 'Y'      % Fix offsets between files when reading in based on 5 min ends of data

check_data_percent ='N' % Check to see how much data exists in file; move to error folder and ignore if below threshold
min_threshold_good=0.98 % keep file if >99% good data

remove_outliers='Y'     % replace outliers with nans? This uses moving window to identify outliers.... makes it much slower.
outlier_stds=10          % use this number of standard deviations from mean to identify bad data [10 is good...needs to be big to avoid throwing out real waves]
out_window = [600 600]  % moving window definition [behind infront] [600 600] is centred 10-minute window with 2 Hz data. If only single number, window is balanced
replace_nans='Y'        % replace nans  by interpolation
maxgap_secs=3700        % maximum gap in seconds to interpolate across

plot_combined = 'Y'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      READ DATA FILES (loop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir(output_directory)
mkdir(error_folder)

%% read chari annual data (strling)
% 
% % filename='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2020_Stirling/StirlingData2020.csv'
% % 
% % t=readtable(filename);
% % ts=datenum(t.Timestamp_AWST_);
% % datestr(ts,'yyyy-mm-dd HH:MM:SS.FFF');
% % 
% % figure; 
% % plot(ts+12/24,t.UnderwaterPressure_m_H20_-10,'color',[.7 .7 .7]);
% % datetick;
% % hold on
% 
% 
% annual_data='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2020_Stirling/StirlingData2020.mat'
% 
% load(annual_data);
% 
% time_adjustment_hours= 12 %19; % addded this value to make data in correct local time
% atime=datenum(t.Timestamp_AWST_(681:end))+time_adjustment_hours/24; % 19 hours error correction don't know why
% ap=t.UnderwaterPressure_m_H20_(681:end);
% 
% [ay am ad ah ami as]=datevec(atime);


%% read annual csv (parmelia)

 filename='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2020_Parmelia/ParmeliaData2020.csv'
% 
t=readtable(filename);
% save('/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2020_Parmelia/ParmeliaData2020.mat','t','filename');
annual_data='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2020_Parmelia/ParmeliaData2020.mat'
% load(annual_data)
ts=datenum(t.Timestamp_AWST_);
% datestr(ts,'yyyy-mm-dd HH:MM:SS.FFF');

figure; 
plot(ts+12/24,t.UnderwaterPressure_m_H20_-10,'color',[.7 .7 .7]);
datetick;
hold on
load fmrawtg2020  %local time % time Sealevel
plot(time,Sealevel+2.4,'k','linewidth',2)

time_adjustment_hours = 12 %19; % addded this value to make data in correct local time
atime=datenum(t.Timestamp_AWST_(1395:end))+time_adjustment_hours/24; % +12 hours error correction don't know why
ap=t.UnderwaterPressure_m_H20_(1395:end);

% do extra adjustment on pre-4 august ( data need +16 hours before 4 August gap)
extra_time_adj_cutoff=datenum(2020,8,4,10,0,0);
id=find(atime<extra_time_adj_cutoff);
% ttmp=atime(id)+4/24; % add extra 4 hours
% ptmp=ap(id);
% figure; 
% hold on; 
% plot(ttmp,ptmp-10,'color',[.7 .6 .6])
% plot(time,Sealevel+2.4,'k','linewidth',2)
extra_time_adjustment_hours = 4
atime(id)=atime(id)+extra_time_adjustment_hours/24; % add extra 4 hours

% keep only the good data

figure; 
hold on; plot(atime,ap-10,'color',[.7 .9 .9])
plot(time,Sealevel+2.4,'k','linewidth',2)
datetick

[ay am ad ah ami as]=datevec(atime);

clear t ts
%% read and clean monthly chunks of data


co=0;

for year = years2do
    d=[];
    tic
    for month = months2do
        co=co+1;
        
        
        %         get indices for current month
        
        mi=find(ay==year & am==month);
        if isempty(mi)
            disp(['No data for ' num2str(year) '/' num2str(month)])
            continue
        end
        
        
        dat=[];ts=[];p=[];
        
        %create a logfile for each month
        
        fid = fopen([output_directory filesep num2str(year) '-' sprintf('%02d',month) '_processing_log.txt'], 'a');
        if fid == -1
            error('Cannot open log file.');
        end
        yourMsg = ' Starting processing.'
        fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
        % fclose(fid);
        
        
        %%-------------------------------------------------------------------------%
        %                       read data
        %%-------------------------------------------------------------------------%
        
        %save logfilename [before loop breaks or removes error files]
        logfilename{co} = annual_data(end-19:end);
        
        
        clear ts p
        
%         subset for the month
        ts=atime(mi);
        p=ap(mi);
        
        
        % ------------------------------------------------------------------------ %
        
        
        
        switch remove_outliers
            case 'Y'
                %flag bad data with NaNs
                %                     badi=find(p<(mean(p,'omitnan')-3*std(p,'omitnan')) | p>(mean(p,'omitnan')+outlier_stds*std(p,'omitnan'))); % to clean up data (simple way but no good)
                [badi] = get_outliers_yh(p,outlier_stds,out_window);
                
                % try to get outliers when instrument changed over
                [badi2]=p<median(p,'omitnan')-2;
                
                if ~isempty(badi)
                    p(badi)=NaN;
                    p(badi2)=NaN;
                    %                         disp([num2str(sum(badi)+sum(badi2)) ' outliers > ' num2str(outlier_stds) ' standard deviations from mean replaced with NaNs']);
                    %
                    yourMsg = [datestr(mean(ts,'omitnan'),'yyyy-mm') ' --  ' num2str(sum(badi)+sum(badi2)) ' outliers > ' num2str(outlier_stds) ' standard deviations from mean replaced with NaNs'];
                    warning(yourMsg)
                    fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                end
        end
        
        
        
%         rename
        combined_ts = ts;
        combined_p = p;
        
        
        
        
        %%
        
        % get unique
        [C,IA,IC] = unique(combined_ts);
        
        % sort
        combined_time= combined_ts(IA);
        combined_p= combined_p(IA);
        
        % get rid of strange situation where time ==0
        badi=find(combined_time==0);
        if ~isempty(badi)
            disp('Removing times = 0');
            combined_time(badi)=[];
            combined_p(badi)=[];
        end
        
%         ADD IN CHECK FOR TIMES = NAN
        
        
        
        figure; plot( combined_time, combined_p);set(gca,'fontsize',14); datetick('x','mm/dd HH:MM','keeplimits')
        
        data=[combined_time combined_p];
        
        % save
        
        if strcmpi(replace_nans,'Y') | strcmpi(replace_nans,'remove_outliers','Y') | strcmpi(fix_offset,'Y')
            processing_applied='Cleaned';
        else
            processing_applied='Raw';
        end
        
        README_time_adjustment=['Time was corrected to local AWST  by adding ' num2str(time_adjustment_hours) ' hours; however prior to '  datestr(extra_time_adj_cutoff) ' an extra ' num2str(extra_time_adjustment_hours) ' hours needed to be added making it +' num2str(time_adjustment_hours+extra_time_adjustment_hours) ' hours total offset'];
        outfname=[output_directory filesep num2str(year) sprintf('%02d',month) '_' location '_' processing_applied '_pressure_data'];
        save(outfname,'data','sRateHz','fix_offset','check_data_percent','min_threshold_good','remove_outliers','outlier_stds','out_window','replace_nans','maxgap_secs','README_time_adjustment');
        
        
        
        
        disp(['Pressure data saved to mat file: ' outfname])
        
        yourMsg =['Pressure data saved to mat file: ' outfname];
        fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg)
        fclose(fid);
        toc
        
        switch plot_combined
            case 'Y'
                f=figure;
                f.Position=[261 648 1818 476];
                plot(combined_time,combined_p,'.');
                datetick('x','yyyy-mm-dd','keeplimits');
                set(gca,'fontsize',14)
        end
        
    end
    
    
    
end



