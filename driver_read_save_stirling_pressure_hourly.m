% Yasha Hetzel 9 Mar 2023
% driver_read_save_stirling_pressure_v9_hourly.m

% Read and save raw data for chari in monthly chunks from HOURLY Cockburn Sound Stirling channel
% marker data from Fremantle Ports (prior to 2021)
% input files are hourly log files with 2Hz sample rate
% should run from within yasha-RBR_tools_v8 directory
% Requirements:
% import_pressure_logfiles.m

% Notes
% - Should probably add a test to make sure the timeseries is continuous --
% - Need to add the data cleaning options contained in the daily version

% UPDATES
% 20230414 
% - saved as v9 for consistency with other working versions 
% - updated only header info



clear all;
clc;
close all;



restoredefaultpath % this just fixes if there are any conflicting functions


% % should work if you run it from the yasha_RBR_tools_v9 directory
% these are not required:
% addpath(genpath(pwd))
% addpath(genpath('rbr-rsktools'))

%%-------------------------------------------------------------------------%
%%                          Settings
%%-------------------------------------------------------------------------%

years2do = [2020] % [2017 2018] % can be vector of years (for loop)
months2do = [8] %[7 8]   % can be vector of months (e.g. [1:12] for loop)

input_directory='/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/2020-08-31T14-00'
output_directory=['/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test']

% move problematic files to here:
error_folder='/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/2020-08-31T14-00/error_folder'

sRateHz=2; % sampling rate in Hz

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
        d=dir([fdir filesep num2str(year) '-' sprintf('%02d',month) '*.log'])
        
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
        
        
        for fi=1:length(d);
            
            co=co+1;
            disp([num2str(year) '-' num2str(month) ' Reading ' num2str(co) ' of ' num2str(length(d)) ' files'])
            lastwarn('');% Clear last warning message
            filename=[d(fi).folder filesep d(fi).name]
            
            
            %%-------------------------------------------------------------------------%
            %                       read data
            %%-------------------------------------------------------------------------%
            %             save logfilename [before loop breaks or removes error files]
            logfilename{co} = d(fi).name ;
            
            
            clear ts p
            [dat ts p] = import_pressure_logfiles(filename);
            %             p=p-10.1325; % remove mean atm pressure. It would be better if you had observations. I have them for this period but there are a few gaps so probably not best. If you are just calculating waves it doesnt matter.
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                yourMsg =[d(fi).name ' is either bad or missing...Skipping and moving to error_folder/'];
                warning(yourMsg)
                fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                movefile(filename,[error_folder filesep])
                continue;
                
            end
            
            
            % check for NaNs
            percgood=sum(~isnan(p))./length(p);
            warning on
            if percgood<1.0 & percgood>0.99
                yourMsg =[d(fi).name '--Only ' sprintf('%0.2f',100*percgood) ' % good data...replacing NaNs with mean'];
                warning(yourMsg)
                fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
                p(isnan(p))=mean(p,'omitnan');
            elseif percgood<0.99
                yourMsg =[d(fi).name '--Only ' sprintf('%0.2f',100*percgood) ' % good data in: ' d(fi).name ' Skipping this file!!!!'];
                warning(yourMsg)
                fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg)
                continue
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
        
        raw_data=[combined_time combined_p];
            
        % save
        outfname=[output_directory filesep num2str(year) sprintf('%02d',month) '_raw_pressure_data'];
        save(outfname,'raw_data');
        
        
        disp(['Raw pressure data saved to mat file: ' outfname])
        
        yourMsg =['Raw pressure data saved to mat file: ' outfname];
        fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg)
        fclose(fid);
        toc
        
   
        
    end
end

