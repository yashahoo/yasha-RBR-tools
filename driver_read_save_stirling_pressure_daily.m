% Yasha Hetzel 13 April 2023
% driver_read_save_stirling_pressure_v9_daily.m

% Read and save raw DAILY [2021+] text files data for chari in monthly chunks from
% Cockburn Sound Stirling channel marker data from Fremantle Ports
% input files are DAILY .txt files with 2Hz sample rate. Contains options
% to clean the data and remove outliers. fix offsets
% Should run from within yasha-RBR_tools_v9 directory
% Requirements:
%
% -import_daily_pressure.m
% -fillgapYLH.m
% -interp_r.m
% -get_outliers_yh.m

% Notes-
% Should probably add a test to make sure the timeseries is continuous; and
% to make perfect time vector
%
% UPDATES
%
% 20230414 yh
% - saved as v9 for consistency with other working versions
% 20230425 yh
% - added cleaning and checking for NaN times
% 20230619 yh
% - added abiity to check for long time gaps for each file and then adjust
% using long term mean rather than 5 mins


clear all;
clc;
close all;



restoredefaultpath % this just fixes if there are any conflicting functions


% % should work if you run it from the yasha_RBR_tools_v9 directory
% these are not required:
addpath(genpath(pwd))
% addpath(genpath('rbr-rsktools'))

%%-------------------------------------------------------------------------%
%%                          Settings
%%-------------------------------------------------------------------------%

years2do = [2022] % [2017 2018] % can be vector of years (for loop)
months2do = [6] %[7 8]   % can be vector of months (e.g. [1:12] for loop)


site= 'Success' %'Stirling % this is important for reading files and naming

% input_directory='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2022_Stirling'
% output_directory=['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2022_Stirling_CLEANED_v2']
% % output_directory=['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2022_Stirling_RAW']
% 
% input_directory='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2021_Stirling'
% output_directory=['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2021_Stirling_CLEANED_v1']

input_directory='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/20230606_WAMSI_FPA_Success_Raw_Digiquartz_Data_2021-23/CombinedDataFeed'
output_directory=['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/Data_2021-2023_Success_CLEANED_v1']


% move problematic files to here:
% % error_folder='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2022_Stirling/error_folder'
% error_folder='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/Data_2021_Stirling/error_folder'
error_folder='/Volumes/Margs_Clone/Pressure_sensor_data/Fremantle_Ports/20230606_WAMSI_FPA_Success_Raw_Digiquartz_Data_2021-23/CombinedDataFeed/error_folder'

sRateHz=2; % sampling rate in Hz

fix_offset = 'Y'      % Fix vertical offsets between files when reading in based on 5 min ends of data or long term mean

check_data_percent ='Y' % Check to see how much data exists in file; move to error folder and ignore if below threshold
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

for year = years2do
    d=[];
    tic
    for month = months2do
        
        
        
        fdir=input_directory; %rename
        %         d=dir([fdir filesep num2str(year) '-' sprintf('%02d',month) '*.log']) % this for hourly log files
        
        d=dir([fdir filesep num2str(year) sprintf('%02d',month) '*.txt']) % for daily txt files
% %         d=d(14:20) %yhtest
        if isempty(d)
            disp(['No files for: ' num2str(year) '-' sprintf('%02d',month) ])
            continue
        end
        
        
        
        dat=[];ts=[];p=[];co=0;
        
        %create a logfile for each month
        
        fid = fopen([output_directory filesep num2str(year) '-' sprintf('%02d',month) '_processing_log.txt'], 'a');
        if fid == -1
            error('Cannot open log file.');
        end
        yourMsg = ' Starting processing.'
        fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
        % fclose(fid);
        
        
        % make data matrices to fill
        %         logfilename= cell(length(d),1); % this doesnt work!!
        logfilename= repmat({''},length(d),1);
        number_samples= NaN(length(d),1);
        num_infiles=length(d);
        
        %         5356800
        combined_ts = NaN(length(num_infiles+1)*24*3600*sRateHz,1);
        combined_p = NaN(length(num_infiles+1)*24*3600*sRateHz,1);
        
%%       
        for fi=1:length(d);
            
            co=co+1;
            disp([num2str(year) '-' num2str(month) ' Reading ' num2str(co) ' of ' num2str(length(d)) ' files'])
            lastwarn('');% Clear last warning message
            filename=[d(fi).folder filesep d(fi).name]
            
% % % %             if fi==14 %test
% % % %                 keyboard
% % % %             end
% % % %             
            %%-------------------------------------------------------------------------%
            %                       read data
            %%-------------------------------------------------------------------------%
            
            % keep last 5 min mean in order to see if there is an offset between files
            endi=[];
            endi=find(~isnan(p),1,'last');
          
            if ~isempty(endi)
                last5mean=mean(p(endi-sRateHz*300:endi),'omitnan');
                last5ts=ts(endi);
                daymean=mean(p,'omitnan');
                longmean=mean(combined_p(combined_p~=0),'omitnan'); % need to account for zeros in data
            else
                last5mean=NaN;
                last5ts=NaN; % just set to NaN
                daymean=mean(p,'omitnan');
                longmean=mean(combined_p(combined_p~=0),'omitnan'); % need to account for zeros in data
            end
            
            
            %             save logfilename [before loop breaks or removes error files]
            logfilename{co} = d(fi).name ;
            
            
            clear ts p
            %             [dat ts p] = import_pressure_logfiles(filename); % this for hourly logfiles
            %             p=p-10.1325; % remove mean atm pressure. It would be better if you had observations. I have them for this period but there are a few gaps so probably not best. If you are just calculating waves it doesnt matter.
            
            switch site
                case 'Stirling'
            [dat ts p] = import_daily_pressure(filename); %
                case 'Success'
            [dat ts p] = import_daily_pressure_success(filename); %
            end
            
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                yourMsg =[d(fi).name ' is either bad or missing...Skipping and moving to error_folder/'];
                warning(yourMsg)
                fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                movefile(filename,[error_folder filesep])
                continue;
                
            end
            % ------------------------------------------------------------------------ %
            %                     fix problem with offset
            % ------------------------------------------------------------------------ %
            
            switch fix_offset
                case 'Y'
                    
                    % catch case where this is first file to define if offset
                    
                    sti=find(~isnan(p),1,'first');
                    first5mean = mean(p(sti:sRateHz*300),'omitnan');
                    
                    % if big gap
                    if ~isnan(last5ts);
                    gapts=ts(sti)-last5ts;
                    is_big_gap=gapts>2/24; % gap is bigger than 2 hours
                    else
                        is_big_gap=logical(1); % say there is a big gap if its teh first one or last time is nan 0?
                    end
                    
                    if isnan(last5mean) & ~is_big_gap % if no big gap and its first data, use it.
                        last5mean = first5mean;

                        yourMsg =[d(fi).name ' - last5mean not available, assuming no offset and getting mean of FIRST 5 mins data'];
                        warning(yourMsg)
                        fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                        
                    elseif  is_big_gap & ~isnan(longmean)  % if there is big gap and you can use long mean
                            last5mean = longmean;
                        yourMsg =[d(fi).name ' - last5mean available, BUT there is big gap (' num2str(gapts*24) ' hrs). Using LONGMEAN (' num2str(longmean) ') for offset'];
                        warning(yourMsg)
                        fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                            
                    elseif ~isnan(last5mean) & is_big_gap & isnan(longmean) % if there is big gap and its the first data, just use beginning height
                            last5mean = first5mean;
                        yourMsg =[d(fi).name ' - last5mean available, BUT there is big gap (' num2str(gapts*24) ' hrs). NO LONGMEAN available...assuming no offset and getting mean of FIRST 5 mins data'];
                        warning(yourMsg)
                        fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                        
                    end
                    
                    % find if big offset between files
                    foffset=first5mean-last5mean;
                    if abs(foffset)>0.3 % changed from 0.5
                        last5mean=last5mean-foffset;
                        yourMsg =[d(fi).name ' - subtracting  ' num2str(foffset) ' to match with previous file'];
                        warning(yourMsg)
                        fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                        
                        p=p-foffset;
                        
                    end
                    
            end
            % ------------------------------------------------------------------------ %
            
            
            
            switch remove_outliers
                case 'Y'
                    %flag bad data with NaNs
                    %                     badi=find(p<(mean(p,'omitnan')-3*std(p,'omitnan')) | p>(mean(p,'omitnan')+outlier_stds*std(p,'omitnan'))); % to clean up data (simple way but no good)
                    [badi] = get_outliers_yh(p,outlier_stds,out_window);
                    
                    %                     try to get outliers when instrument changed over
                    [badi2]=p<median(p,'omitnan')-2;
                    
                    if ~isempty(badi)
                        p(badi)=NaN;
                        p(badi2)=NaN;
%                         disp([num2str(sum(badi)+sum(badi2)) ' outliers > ' num2str(outlier_stds) ' standard deviations from mean replaced with NaNs']);
%                         
                          yourMsg = [d(fi).name '--  ' num2str(sum(badi)+sum(badi2)) ' outliers > ' num2str(outlier_stds) ' standard deviations from mean replaced with NaNs'];
                                warning(yourMsg)
                                fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                    end
            end
            
            switch check_data_percent
                case 'Y'
                    
                    % check for NaNs (after filling small gaps)
                    percgood=sum(~isnan(p))./length(p);
                    warning on
                    if percgood<1.0 & percgood>min_threshold_good % percgood<1.0 & percgood>0.99
                        yourMsg =[d(fi).name ' -- ' sprintf('%0.4f',100*percgood) ' % good data'];
                        disp(yourMsg)
                        fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                        %                         p(isnan(p))=mean(p,'omitnan');
                        switch replace_nans
                            case 'Y'
                                % interpolate across small gaps only
                                maxgap=maxgap_secs*sRateHz;
                                [pfilled,gaps_filled,gaps_unfilled,bind,gapindxj0,gapindxj1]=fillgapYLH(p,maxgap);
                                disp([num2str(gaps_filled) ' gaps filled and ' num2str(gaps_unfilled) ' gaps NOT filled'])
                                p=pfilled;
                                
                                yourMsg = [d(fi).name '--  ' num2str(gaps_filled) ' gaps (< ' num2str(maxgap) ') filled and ' num2str(gaps_unfilled) ' gaps NOT filled'];
                                warning(yourMsg)
                                fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                                
                        end
                        
                    elseif percgood<min_threshold_good
                        yourMsg =[d(fi).name '--Only ' sprintf('%0.4f',100*percgood) ' % good data in: ' d(fi).name ' Skipping this file!!!!'];
                        warning(yourMsg)
                        fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg)
                        continue
                    end
            end
            
            
            
            
            
            
            % concatinate the raw data
            if co==1
                starti=1;
            else
                
                starti=starti+length(ts);
            end
            
            % %             fix problem where
            %             if ts(1)==0;
            %                 ts(1)=[];
            %                 p(1)=[];
            %             end
            
            try
                combined_ts(starti:starti+length(ts)-1) = ts;
                combined_p(starti:starti+length(ts)-1) = p;
                
            catch
                keyboard
                yourMsg =[d(fi).name ' Failed to concatinate raw data to monthly  ...skipping and moving to error_folder/'];
                warning(yourMsg)
                fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                movefile(filename,[error_folder filesep])
                continue;
                
            end
            
        end   % loop through files
        

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
        
        %         also remove where NaNs at end
        badit=find(isnan(combined_time));
        if ~isempty(badit)
            disp('Removing times = NaN');
            combined_time(badit)=[];
            combined_p(badit)=[];
        end
        
        
        figure; plot( combined_time, combined_p);set(gca,'fontsize',14); datetick('x','mm/dd HH:MM','keeplimits')
        
        data=[combined_time combined_p];
        
        % save
        
        if strcmpi(replace_nans,'Y') | strcmpi(replace_nans,'remove_outliers','Y') | strcmpi(fix_offset,'Y')
           processing_applied='Cleaned';
        else
            processing_applied='Raw';    
        end
        
        outfname=[output_directory filesep num2str(year) sprintf('%02d',month) '_' processing_applied '_pressure_data'];
        save(outfname,'data','sRateHz','fix_offset','check_data_percent','min_threshold_good','remove_outliers','outlier_stds','out_window','replace_nans','maxgap_secs');
        
 
        
        
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

