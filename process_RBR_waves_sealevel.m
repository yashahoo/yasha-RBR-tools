% Yasha Hetzel April 2021
% process_RBR_waves_sealevel_v9.m
% This m-file reads raw .rsk file downloaded from RBR pressure sensors
% using rsk-tools (downloadad from rbr website)
% and saves to netcdf file. then it reads the netcdf file and calculates
% waves (optional), tidal analysis, lowpass filter, highpass filter (optional)
% and makes a number of plots that are saved to a directory along
% with the processed data as matfile (optional) or the hourly wave data as a csv text file
% **There are a lot of m-files required so you need to set the path to make
% sure all the m-files in the yasha_RBR_tools are available.
% 
%  Waves now correctly calculated with pressure response for depth using
%  Oceanlyz toolbox
%  http://zarmesh.com/wp-content/uploads/2016/05/Wind-wave-analysis-in-depth-limited-water-using-OCEANLYZ-A-MATLAB-toolbox.pdf
%  https://oceanlyz.readthedocs.io/en/latest/
% 
% UPDATES
%  20220913 yh
% fixed saving to not require export_fig
% added subsample averaging interval < 1 hr plottting and csv export
% changed number samples for IG waves processing to be 4096 not 2048
% fixed deleting of existing nc files so doesn't cause problems
% noted that location string not used in file naming (add to TODO list)

% 20220920 yh
% added opetion to correct pressure for depth
% added write trimmed data to nc file
% added write pressure corrected data to nc file
% added offset from UTC in read_RBRyh.m
% added offset from UTC
% 20220928 yh
% added option to  adjust to Chart Datum

% 20221019 yh
%  - this version uses oceanlyz function to correct pressure response and to calculate waves
%  - option to use zero-crossing or spectral methods and to separate into sea and swell if using spectral method
%  - read_RBRyh.m has been adjusted to not remove mean atmospheric pressure
%  - uses bom observations of mslp to remove atmosphere if available
%  - retains ability to choose chart datum (if you have correct measurements) or mean sealevel datum

% 20230412 yh
% -modified to use read_RBRyh_v2.m which correctly reads pressure and temperature channels in new RBR sensor

% 20230414 yh
% - saved as v9 for consistency with other working versions
% - v9 and process_RBR_waves_sealevel are the same

%%


clear all;
clc;
close all;

% required [contained in yasha_RBR_tools]:
% t-tide
% rotateXlabels
% subaxis
% exportfig
% lanczosfilter
% oceanlyz toolbox [optional] https://oceanlyz.readthedocs.io/en/latest/
% and some others including rbr rsk-tools


restoredefaultpath % this just fixes if there are any conflicting functions

% set the paths here relevant to your computer
% addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/yashatools/yasha_RBR_tools')); % [edit]  this should contain all the functions you need, except the rsk-tooks that are required to read the raw data (see below)
% % addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rbr-rsktools-2.3.0')); %  [edit] this should be installed on your computer somewhere (it has been tested ususing old version)
% addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rbr-rsktools')); % use this old version to ensure it reads temperature!

% this should work if you run it from the yasha_RBR_tools directory
addpath(genpath(pwd))
% addpath(genpath('rbr-rsktools'))
%%-------------------------------------------------------------------------%
%% -------------------------settings-------------------------------------- %
%%-------------------------------------------------------------------------%

% 

 
%%-------------------------------------------------------------------------%
% set location and name of raw rsk file. by default, figs will be saved to same directory
%%-------------------------------------------------------------------------%
% location='Port_Geographe'; rawfile='/Users/00068592/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/SN206857_Dirk_Port_Geographe_20210831-20211021/206857_20220629_0954.rsk'

% % Deployment 1 - StuBarr
% datum_correction=[1.54-0.95] ;site_depth=1.5;  location='Abbey'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/20220818_Busselton_RBR_ABB_DUNS_EB_StuBarr/ABBEY_206860_20220818_0948.rsk'
% datum_correction= [1.74-1.2]  ;site_depth=1.5;  location='OldDunsborough'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/20220818_Busselton_RBR_ABB_DUNS_EB_StuBarr/DUNS_206857_20220818_1312.rsk' %done
% datum_correction= [1.9-1.5] ; site_depth=1.5;  location='EastBusselton'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/20220818_Busselton_RBR_ABB_DUNS_EB_StuBarr/EASTBUSSO_209847_20220818_1102.rsk' % done

% site_depth=6.5; location='BusseltonJetty341_end'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/208787_20220911_1406.rsk' % end of jetty
% site_depth=1.0; location='BusseltonJetty1'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/077826_20220911_1308.rsk' %location unknown
% site_depth=5.5;  location='BusseltonJetty204'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/208788_20220911_1353.rsk' % middle
% site_depth=3.0;  location='BusseltonJetty147A'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/208790_20220911_1444.rsk' % outside sandbar
% site_depth=3.0;  location='BusseltonJetty133A';  rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/209845_20220911_1418.rsk' % centre sandbar
% site_depth=3.0; location='BusseltonJetty108A'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/124227_20220911_1459.rsk' % 108A inside sandbar
% site_depth=3.0; location='BusseltonJetty26A'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/124027_20220911_1321.rsk' % 26A nearshore ( this appears to be corrupt or old firmware. must check before re-deploying if update possible). CANNOT READ IT.

% Deployment 2 - StuBarr
% datum_correction=[1.54-0.95] ; site_depth=1.5;  location='Abbey'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/20221027_Busselton_RBR_ABB_DUNS_EB_StuBarr/ABBEY_D2_206860_20221027_1035.rsk'
% datum_correction= [1.74-1.2]  ; site_depth=1.5;  location='OldDunsborough'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/20221027_Busselton_RBR_ABB_DUNS_EB_StuBarr/DUNS_D2_206857_20221027_1254.rsk' %done
% datum_correction= [1.9-1.5] ;     site_depth=1.5;  location='EastBusselton'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/20221027_Busselton_RBR_ABB_DUNS_EB_StuBarr/EastBusselton_D2_209847_20221027_1053.rsk' % done

% Deployment 5 - StuBarr
% datum_correction=[1.54-0.95] ; site_depth=1.5;  location='Abbey'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/SN_124229_StuBarr_Abbey_Busselton_20210206-20210504/ABBEY_D5_124229_20230619_0854.rsk'

% datum_correction= [NaN] ; site_depth=1.5;  location='Example'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/20220818_Busselton_RBR_ABB_DUNS_EB_StuBarr/EASTBUSSO_209847_20220818_1102.rsk' % done

% datum_correction= [NaN] ; site_depth=10;  location='Stirling'; rawfile='/Volumes/Margs_Clone/RBRpressure_sensors/SN_206854_20230403_chari_new/206854_20230403_1218.rsk' % done 6 depth


% % datum_correction= [NaN] ; site_depth=7;  location='SB1'; rawfile='/Volumes/YASHA16/124227_20230523_1312.rsk'
% datum_correction= [NaN] ; site_depth=7;  location='SB2'; rawfile='/Volumes/YASHA16/208790_20230523_1259.rsk'


datum_correction=[NaN] ; site_depth=2.8;  location='Gracetown_wavebuoy'; rawfile='/Volumes/Margs_Clone/Pressure_sensor_data/RBRpressure_sensors/SN124229_Gracetown_Wavebuoy_mooring_20230827/Gracetown_Wavebuoy_124229_20230828_1021.rsk'

%%-------------------------------------------------------------------------%

% EXAMPLE measurements required to change datum from MSL to Chart datum
% Dunsborough:
% start_WL_datum_ht=[1.74-1.2+0.6] ;% [NaN] % [ht of deck - rbr to deck + rbr to water surface]
% Deck Elevation - 1.74mCD
% Distance between RBR and Deck - 1.2m
% D1 RBR to WL - 0.6m 
% D2 RBR to WL - 0.6m 
% Port Geographe Tide Gauge - 1.1mCD
% Time: 13:40
%
% Abbey:
% start_WL_datum_ht=[1.54-0.95+0.66] ;% [NaN] % [ht of deck - rbr to deck + rbr to water surface]
% Deck - 1.54mCD
% Distance between RBR and Deck - 0.95m
% RBR to WL - 0.66m
% D2 RBR to WL - 0.28m 
% Port Geographe Tide Gauge - 1.25mCD
% Time: 12:30
%
% East Busselton:
% start_WL_datum_ht=[1.9-1.5+1]
% Deck Elevation - 1.9mCD
% Distance between RBR and Deck - 1.5m
% RBR to WL - 1m
% D2 RBR to WL - 0.32m 
% Port Geographe Tide Gauge - 1.4mCD
% Time: 11.00am
%%-------------------------------------------------------------------------%
% % deployment info
%%-------------------------------------------------------------------------%
sample_rate_hz=2; % sampling in hz (normally 2)
subsample_mins=1; % water level averaging interval in minutes (all other calculations done in hourly intervals)

% provide absolute levels if you want to convert to Chart Datum:
adjust_to_chart_datum = 'N'             %
% datum_correction= [1.74-1.2]; % [ht of deck - deck to RBR]
% % % % [no longer used: start_WL_datum_ht=[1.54-0.95+0.66] ; %[1.54-0.95+0.66] %[1.74-1.2+0.6] %[1.54-0.95+0.66] ;% [NaN] % [ht of deck - rbr to deck + rbr to water surface] %round_start_finish must  be = 'N' !!!;  this to be added to the MSL zero height for first 30*sample_rate_hz observations. Make Nan if unknown or don't want to be chart datum

%%-------------------------------------------------------------------------%
% OPTIONS -   'Y' if you want to calculate/plot; 'N' if don't
%%-------------------------------------------------------------------------%
round_start_finish = 'N'            % if you want nice start and end times (e.g. on the 1/2 hour). Must not round if adjusting to CHart Datum!!!
calculate_waves    = 'Y'            % calculate wave parameters including Hs, Tz, etc.   
calculate_IG       = 'Y'            % calculate Infragravity waves (shoudl be 'Y' if calculating waves)
get_time_freq      = 'Y'            % use chari's function to calculate time-frequency plot
do_highpass_filter = 'N'            % do highpass filter to remove waves > max_period
plot_highpass_weeks= 'N'            % plot weekly highpass filter (was useful to look at boatwakes in albany)
save_matfile       = 'N'            % save the processed data to matfile
save_matfig        = 'N';           % save the .fig files?
print_fig          = 'Y';           % print to image file?
use_export_fig     = 'Y'            % if print_fig ='Y' -->  'N' uses native matlab print function (default) OR 'Y' uses export_fig to make nice plots
img_type           = {'pdf','png'}; % if print_fig ='Y' -->  {'pdf','png','eps'} % options. can be multiple. pdf will not save if native matlab unless edit function

%%-------------------------------------------------------------------------%
%  Use observations of air pressure to correct data? [recommended if possible]
%%-------------------------------------------------------------------------%
remove_dynamic_atm = 'N'            % Y if you have mslp observations to correct data... if N it will use subtracted 10.1325 db
% load air pressure data [ and get mslp data for closest station)
% must manualy define which station to use (e.g. 26 below is Busselton Airport)
 bom_index=26% % 26;
% 1. Witchcliffe West
% 2. Eucla
% 3. Red Rocks Point
% 4. Rowley Shoals
% 5. Adele Island
% 6. Browse Island
% 7. Broome
% 8. Port Hedland
% 9. Karratha
% 10. Bedout Island
% 11. Learmonth
% 12. Onslow Airport
% 13. Barrow Island
% 14. Carnarvon
% 15. Shark Bay Airport
% 16. North Island
% 17. Geraldton Airport
% 18. Perth Airport
% 19. Rottnest Island
% 20. Ocean Reef
% 21. Swanbourne
% 22. Garden Island
% 23. Hillarys Point Boat Harbour
% 24. Cape Leeuwin
% 25. Cape Naturaliste
% 26. Busselton Airport
% 27. Esperance
% 28. Hopetoun North
% 29. Bunbury
% 30. Mandurah
% 31. Walpole North
% 32. Albany Airport

%%-------------------------------------------------------------------------%
%      Wave calculation settings [for oceanlyz toolbox]
%               if calculate_waves = 'Y' above
%             default settings can be left 'as is'
%%-------------------------------------------------------------------------%
InputType='pressure' %'waterlevel' %'pressure'   %'waterlevel' %           % Data input ['pressure'] or  'waterlevel'. If pressure, the signal attenuation with depth will be applied [recommended].
OutputType='wave+waterlevel' %'wave' %'wave+waterlevel'      % ['wave'], ['wave+waterlevel'];
AnalysisMethod=landcruiser co'zerocross' %'spectral' %'zerocross'  %'zerocross' %'spectral'         % Wave calculation method. 'zerocross' or ['spectral']. use zerocross if you want Hs, but very similar to Hm0 and swell/sea are useful
burst_duration=3500;              % Interval of time window (in seconds) over which to calculate waves
fs=sample_rate_hz;                % Frequency of data in Hz
heightfrombed= 0 %site_depth-1 ;  % [Default=0]   % height of instrument above bed (m). Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'.If unknown assume = 0
dispout='no';                     % Display output of calculations at each burst ('yes or 'no')
Rho=1024;                         % Seawater density (Varies)
SeparateSeaSwell='yes'            % 'yes or 'no' (only works with spectral)
fmaxswell=1/8 %0.25;              % Maximum swell frequency (spectral only). 1/Period (ie. minimum period a swell can have). It should be between 0 and (fs/2). Only required if SeparateSeaSwell='yes' and AnalysisMethod='spectral'
fpminswell=1/25 %0.1;             % Minimum swell frequency (spectral only)
max_period= 25;                   % max period to use in wave calculations and filtering
%-------------------------------------------------------------------------%
%   [optional] set time ranges or flag bad data range
%-------------------------------------------------------------------------%
% time range to process.
% If unknown, make NaN and will be selected interactively [Default]
% begin=datenum(2022,8,18,10,12,32);
% finish=datenum(2022,9,21,12,17,22);
% begin=datenum(2022,7,1,10,54,34);
% finish=datenum(2022,8,16,9,11,31);
begin=NaN;
finish=NaN;

% manually flag bad data with nans [if no bad data make badstart/end = NaN]. this
% seems to cause some problemss,so avoid for now if possible
% badstart=datenum(2018,2,4,4,17,3);
% badend=datenum(2018,2,4,4,19,20);
badstart=NaN;
badend=NaN;
%-------------------------------------------------------------------------%
%        set where to save figures and data
%-------------------------------------------------------------------------%
font_size=12; % set font size for figs
switch calculate_waves
    case 'Y'
whereput=[rawfile(1:end-4) '_' location '_processed_' AnalysisMethod]; % this saves to same directory where rsk file located
    case 'N'
whereput=[rawfile(1:end-4) '_' location '_processed_NO_waves']; % this saves to same directory where rsk file located        
end
mkdir(whereput);

%%-------------------------------------------------------------------------%
%%-------------------------------------------------------------------------%
%%                     get BOM MSLP and wind data
%%-------------------------------------------------------------------------%
%%-------------------------------------------------------------------------%
switch remove_dynamic_atm
    case 'Y'
        disp('Loading Bom mslp data...')
        if ~exist('archive','var')
            try
                load('BOM_obs_archived.mat')
            catch
                disp('Local BOM data does not exist...downloading from yasha website...')
                [WEBDATA]=websave('BOM_obs_archived.mat','https://sealevelx.ems.uwa.edu.au/plot/yasha/BOM_obs_archived.mat')
                load(WEBDATA)
            end
        end
        
%         bom_index=26;
        
        mslp=archive(bom_index).press_msl;
        mslp_time=archive(bom_index).mtime_UTC;
        bom_sta= archive(bom_index).name{1}
end


%%-------------------------------------------------------------------------%
%%-------------------------------------------------------------------------%
%%      read rsk, create netcdf and begin processing
%%-------------------------------------------------------------------------%
%%-------------------------------------------------------------------------%
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
ncref=[whereput filesep filebase '_' location '.nc'];



% check that nc file doesnt yet exist, delete if it does
if isfile(ncref) %exist(ncref,'file') %exist('ncref','var')
    disp('Netcdf output file already exists, delete it?')
    
    promptMessage = sprintf('OVERWRITE netcdf output file or QUIT?');
    promptTitle = sprintf('Netcdf file exists');
    button = questdlg(promptMessage,promptTitle,'Quit','OVERWRITE','Quit');
    if strcmpi(button, 'Quit')
        close all
        return; % Or break or continue
    elseif strcmpi(button, 'OVERWRITE')
        disp(['Deleting: ' ncref])
        delete(ncref)
        clear ncref
    end
    %     disp(['Deleting: ' ncref])
    
end

disp('Previewing rsk file...')
% read raw file and save to netcdf
[ncref]=read_RBRyh_v2(rawfile,location,whereput);
% ncref='/Volumes/Margs_Clone/RBRpressure_sensors/SN202078_Smiths_20210409/202078_20210409_1525_processed/202078_20210409_1525.nc'

%
% % % this used later
% ncref=[whereput '/' rawfile(end-23:end-4) '.nc'];
close all;

% % % % % fix broken time
% % % % load('east_busso_offset.mat')
% % % % ts=ncread(ncref,'time')+datenum(2000,1,1)+dt;
% % % % raw_ts=ncread(ncref,'time')+datenum(2000,1,1)+dt;

ts=ncread(ncref,'time')+datenum(2000,1,1);
raw_ts=ncread(ncref,'time')+datenum(2000,1,1);
raw_press=ncread(ncref,'press');
instrument=ncreadatt(ncref,'/','Serial_Number'); % serial number here for id

try
    Offset_from_UTC=ncreadatt([ncref],'time','Offset_from_UTC')
catch
    disp('Offset_from_UTC not available')
end

try
    % temp=ncread([dir_ref ncref],'temperature');
    temp=ncread([ncref],'temperature');
catch
    disp('Temperature not available')
end


%% correct with mslp observations if available

switch remove_dynamic_atm
    case 'Y'
        %subtract standard atm pressure or use observed mslp data
        if exist('mslp_time','var')
            try
                corrected_with_mslp_obs='Y';
                disp(['Atmospheric pressure has been removed using bom data from: ' bom_sta])
                imslp=interp1(mslp_time,mslp,ts)./100;
                %                 raw_press=p+10.1325; % this was needed for old nc files but now its not added
                %                 raw_press=p;
                p=raw_press-imslp;
            catch
                corrected_with_mslp_obs='N';
                
            end
        else
            corrected_with_mslp_obs='N';
            p=raw_press-10.1325;
            disp('Atmospheric pressure has been removed by subtracting 10.1325 dB')
            
        end
    case 'N'
        corrected_with_mslp_obs='N';
        p=raw_press-10.1325;
        disp('Atmospheric pressure has been removed by subtracting 10.1325 dB')
        
end
%% interactively select start and finish times [OPTIONAL].


% begin=ts(10); finish = nan % yh 20221031
if isnan(begin) | isnan(finish)
    
    promptMessage = sprintf('Do you want to INCLUDE first data point or SELECT both start/finish or use ALL data?');
    promptTitle = sprintf('Selecting start / end');
    button = questdlg(promptMessage,promptTitle,'Include','Select','All','Select');
    if strcmpi(button, 'Include')
        %         close all
        %         return; % Or break or continue
        
        begin=ts(1);
        disp('Starting at index = 1 ');
        if isnan(finish)
            [finish]=get_finish(ts, p);
        end
        
    elseif strcmpi(button, 'Select')
        % interactively select start and finish times
        %         if isnan(begin) && isnan(finish)
        disp('Selecting start and finish');
        [begin finish]=get_beginfinish(ts, p);
        title('Done selecting start + finish')
        %         end
        
        %         disp(['Deleting: ' ncref])
        %         delete(ncref)
        %         clear ncref
    elseif strcmpi(button, 'All')
        
        begin =ts(1);
        finish = ts(end);
        disp('Using  all data');
    end
    
    
end

% old stuff:
% % interactively select start and finish times
% if isnan(begin) && isnan(finish)
% [begin finish]=get_beginfinish(ts, p);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   if you want nice start and end times (e.g. on the 1/2 hour)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if ~isempty(badi)
ts(badi)=[];
p(badi)=[];
else
   disp('Not removing any bad data!!!')
end

if exist('temp','var')
    temp(badi)=[];
end

figure;
hold on
plot(ts,p,'b')
datetick
datetick('x','yyyy-mm-dd HH:MM:SS','keeplimits')

% close

% rename for saving
raw_trimmed_press=p;
trimmed_time=ts;


% set some plot limits
xlims=[ts(1) ts(end)]; % all
xlims=[begin finish];% good
xt=[xlims(1):xlims(2)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% adjust to datum [optional]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start_WL_datum_ht= [ht of deck - rbr to deck + rbr to water surface]
%                    [1.74       -      1.2    +        0.6          ] Duns

switch adjust_to_chart_datum
    case 'Y'
        if ~isnan(datum_correction)
            zero_window_secs=30;
            
            dp=mean(p(1:zero_window_secs*sample_rate_hz))-0;
            
            % %     test=mean(raw_trimmed_press (1:zero_window_secs*sample_rate_hz));
            
%             p_CD=p-dp+start_WL_datum_ht;
%             p=p_CD;
%             datum='CD';
%             datum_correction=-dp+start_WL_datum_ht % add this to raw data, or subtract from  corrected to get raw
            
            p_CD=p+datum_correction;
            p=p_CD;
            datum='CD';
            
            disp('water level /pressure adjusted to chart datum')
            ylabel('Height above Chart Datum (m)')
            
        else
            datum='MSL'
            datum_correction=mean(p, 'omitnan');
        end
        
    case 'N'
        datum='MSL'
        datum_correction=mean(p, 'omitnan');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if bom data available plot wind arrows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch remove_dynamic_atm
    case 'Y'

%         bom_index
        num_days_data=ceil(ts(end)-ts(1));
        
        weeks=floor(ts(1));
        for wi=1:floor(num_days_data/7)-1
        weeks=[weeks; floor(ts(1)+7*wi)];
        end
        
        for k=1:length(weeks)
        [archive]=plot_date_BOM_wind_v2([bom_index],datestr(weeks(k)),7,15,15);
        set(gca,'fontsize',14)
        fname=[whereput '/' instrument '_' location '_WindArrow_week_' num2str(k)];
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
        
        end
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add trimmed and corrected time and pressure to netcdf file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Adding new variables to netcdf file....')

        nccreate(ncref,'corrected_press','Dimensions',{'corrected_time',length(p)});
        ncwrite(ncref,'corrected_press',p);
        ncwriteatt(ncref,'corrected_press','Units','Decibars');
        ncwriteatt(ncref,'corrected_press','Datum',[ datum ' ; correction of ' num2str(datum_correction) ' db added to raw data with atm removed']);

switch remove_dynamic_atm
    case 'N'
        ncwriteatt(ncref,'corrected_press','description','Good pressure data only, mean atmospheric pressure 10.1325 db removed');
    case 'Y'
        ncwriteatt(ncref,'corrected_press','description',['Good pressure data only, and atmospheric mslp removed using used BOM station:  ' bom_sta '.']);
        disp('Adding closest bom observations used to correct data')
        nccreate(ncref , 'mslp','dimensions',{'time',length(imslp)});
        ncwriteatt(ncref,'mslp','units','Millibars');
        ncwriteatt(ncref,'mslp','description',['Mean sea level pressure from closest Bom station: ' bom_sta '.']);
        ncwrite(ncref,'mslp',imslp);
end

        nccreate(ncref , 'corrected_time','Dimensions',{'corrected_time',length(ts)});
        ncwriteatt(ncref,'corrected_time','Units','Days since 2000-01-01T00:00:00');
        ncwriteatt(ncref,'corrected_time','UTC offset',Offset_from_UTC);
        ncwriteatt(ncref,'corrected_time','description','Time corresponding to pressure with atm removed and trimmed to good data. ');
        ncwrite(ncref,'corrected_time',ts-datenum('2000-01-01'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Use OceanLyz wave toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch calculate_waves
    case 'Y'
        % this to be added in next version 7 and chari method removed
        
        n_bursts=floor(length(ts)./(burst_duration*fs));
        % btime=[ts(1):burst_duration/burst_duration/24:ts(n_bursts*burst_duration)];
        % n_bursts=100 %
        % create time vector corresponding to waves
        % burst_times=[ts(1):burst_duration/3600/24:ts(n_bursts*burst_duration*fs-1)];
        % % this is wrong
        
        samples=burst_duration*fs;
        co=0;
        for i=1:n_bursts
            
            i1=co*samples+1;
            i2=i1+samples-1;
            burst_times(i)=mean(ts(i1:i2));
            co=co+1;
        end
        

        

        %Read data
        current_folder=pwd;                  % Current (OCEANLYZ) path
        switch InputType
            case 'pressure'
                water_pressure=p*10000; %              % Load data
            case 'waterlevel'
                water_pressure=p;
        end
        % cd(current_folder)                   % Change current path to OCEANLYZ folder
        water_pressure=water_pressure(1:n_bursts*burst_duration*fs);
        
        
        %add path for OCEANLYZ folder
        addpath(genpath('oceanlyz_2_0'))
        % cd('oceanlyz_2_0')
        
        %Create OCEANLYZ object
        clear ocn %Optional
        ocn=oceanlyz;
        
        
        %Input parameters
        ocn.data= water_pressure;
        ocn.InputType=InputType; %'pressure' %'waterlevel';
        ocn.OutputType=OutputType; %'wave' %'wave+waterlevel';
        ocn.AnalysisMethod=AnalysisMethod; %'spectral';
        ocn.n_burst=n_bursts;
        ocn.burst_duration=burst_duration; %3600;
        ocn.fs=fs;
        %     ocn.fmin=1/20
        %     ocn.fmax=ocn.fs/2;  %0.1250
        ocn.fmin=0.05;                    %Only required if ocn.AnalysisMethod='spectral'
        ocn.fmax=ocn.fs/2;                %Only required if ocn.AnalysisMethod='spectral'
        ocn.fmaxpcorrCalcMethod='auto';   %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
        ocn.Kpafterfmaxpcorr='constant';  %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
        ocn.fminpcorr=0.15;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
        ocn.fmaxpcorr=0.55;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
        ocn.heightfrombed=heightfrombed   % %0.05;           %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
        ocn.dispout=dispout;              %'no'
        ocn.Rho=1024;                     %Seawater density (Varies)
        
        % define inputs for separating sea and swell
        switch AnalysisMethod % cannot separate sea/swell with zero cross method, so check this.
            case 'zerocross'
                SeparateSeaSwell='no';
                ocn.SeparateSeaSwell=SeparateSeaSwell;
            case 'spectral'
                ocn.SeparateSeaSwell=SeparateSeaSwell;
        end
        % if SeparateSeaSwell='yes', define these:
        %maximum swell frequency
        ocn.fmaxswell= fmaxswell ;        %1/8 %0.25;
        %                                 Maximum frequency that swell can have (It is about 0.2 in Gulf of Mexico) in (Hz)
        %                                     It should be between 0 and (fs/2)
        %                                     Only required if SeparateSeaSwell='yes' and AnalysisMethod='spectral'
        
        %Minimum swell frequency
        ocn.fpminswell=fpminswell;       %1/25 0.1;
        %                                 Minimum frequency that swell can have (it is used for Tpswell calculation) in (Hz)
        %                                     It should be between 0 and (fs/2)
        %                                     Only required if SeparateSeaSwell='yes' and AnalysisMethod='spectral'
        
        
        
        
        
        %Run OCEANLYZ
        ocn.runoceanlyz()
        
        % cd ..
% ----------------------------Plots-------------------------------------- %
        % make sure time is correct length
        burst_times=burst_times(1:size(ocn.wave.Burst_Data,1));
        
        
        xlims=[ts(1) ts(end)]; % all
        xlims=[begin finish];% good
        xt=[xlims(1):xlims(2)];
        yt=[0:2:30];
        %
        % sp_corrected1=ocn.wave;
        % sep_swell= ocn.SeparateSeaSwell;
        switch AnalysisMethod
            case 'spectral'
                
                switch SeparateSeaSwell
                    case 'yes'
                        %         figure;
                        %         hold on;
                        %         plot(ocn.wave.Hm0swell);
                        %         plot(ocn.wave.Hm0sea);
                        %         plot(ocn.wave.Hm0,'k');
                        %         legend('swell','sea','Hm0');
                        %
                        %         figure;
                        %         hold on;
                        %
                        %         h1=plot(ocn.wave.Tp);
                        %         h2=plot(ocn.wave.Tpswell);
                        %         h3=plot(ocn.wave.Tpsea);
                        %         legend('Tp','swell','sea');
                        %         ylabel('Period (s)')
                        %
                        %         h1.Color=[0 0 0];
                        %         h2.Color=[ 0    0.4470    0.7410];
                        %         h3.Color=[0.8500    0.3250    0.0980];
                        
                        f=figure;
                        f.Position=[66 207 1600 450];
                        f.Color='w';
                        % set(gcf,'Position',[200 100 700 300])
                        % set(gcf,'PaperPositionMode','auto')
                        plot(burst_times,ocn.wave.Hm0swell,'.') % ,'.'
                        hold on
                        plot(burst_times,ocn.wave.Hm0sea,'.')
                        plot(burst_times,ocn.wave.Hm0,'k')
                        legend('Hm0swell','Hm0sea','Hm0')
                        ylabel('Wave height (m)')
                        set(gca,'xlim',xlims,'xtick',xt)
                        set(gca,'ylim',[0 ceil(max(ocn.wave.Hm0))])
                        datetick('x','mm/dd','keeplimits','keepticks')
                        rotateXLabels( gca(), 90 )
                        xlabel(datestr(mean(burst_times),'YYYY'))
                        title([instrument ' ' location ],'interpreter','none')
                        fname=[whereput '/' instrument '_' location '_WAVE_HEIGHT_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
                        % print('-dpng','-r200', fname);
                        % print('-dpdf', '-painters', fname);
                        %         export_fig(fname,'-pdf','-png','-transparent')
                        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
                        
                        f=figure;
                        f.Position=[200 100 1600 450];
                        f.Color='w';
                        h1=plot(burst_times,ocn.wave.Tp)
                        hold on
                        h2=plot(burst_times,ocn.wave.Tpswell,'.')
                        h3=plot(burst_times,ocn.wave.Tpsea,'.')
                        set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
                        legend('Tp','swell','sea');
                        ylabel('Period (s)')
                        title([instrument ' ' location ],'interpreter','none')
                        datetick('x','mm/dd','keeplimits','keepticks')
                        %         xlabel(datestr(mean(burst_times),'YYYY'))
                        if exist('Offset_from_UTC','var')
                            xlabel([datestr(mean(burst_times),'YYYY') ' (UTC+' Offset_from_UTC ')'])
                        else
                            xlabel(datestr(mean(burst_times),'YYYY'))
                        end
                        rotateXLabels( gca(), 90 )
                        h1.Color=[0 0 0];
                        h2.Color=[ 0    0.4470    0.7410];
                        h3.Color=[0.8500    0.3250    0.0980];
                        
                        fname=[whereput '/' instrument '_' location '_WAVE_PERIOD_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
                        % print('-dpng','-r200', fname);
                        % print('-dpdf', '-painters', fname);
                        %         export_fig(fname,'-pdf','-png','-transparent')
                        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
                        
                        
                    case 'no'
                        
                        %                     figure;
                        %                     hold on;
                        %                     plot(ocn.wave.Hm0);
                        %                     legend('Hm0');
                        %
                        %                     figure;
                        %                     hold on;
                        %                     plot(ocn.wave.Tp);
                        %                     legend('Tp');
                        %                     ylabel('Period (s)')
                        
                        f=figure;
                        f.Position=[66 207 1600 450];
                        f.Color='w';
                        plot(burst_times,ocn.wave.Hm0,'k')
                        legend('Hm0')
                        ylabel('Wave height (m)')
                        set(gca,'xlim',xlims,'xtick',xt)
                        set(gca,'ylim',[0 ceil(max(ocn.wave.Hm0))])
                        datetick('x','mm/dd','keeplimits','keepticks')
                        rotateXLabels( gca(), 90 )
                        xlabel(datestr(mean(burst_times),'YYYY'))
                        title([instrument ' ' location ],'interpreter','none')
                        fname=[whereput '/' instrument '_' location '_WAVE_HEIGHT_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
                        % print('-dpng','-r200', fname);
                        % print('-dpdf', '-painters', fname);
                        %         export_fig(fname,'-pdf','-png','-transparent')
                        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
                        
                        f=figure;
                        f.Position=[200 100 1600 450];
                        f.Color='w';
                        h1=plot(burst_times,ocn.wave.Tp)
                        set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
                        legend('Tp','swell','sea');
                        ylabel('Period (s)')
                        title([instrument ' ' location ],'interpreter','none')
                        datetick('x','mm/dd','keeplimits','keepticks')
                        %         xlabel(datestr(mean(burst_times),'YYYY'))
                        if exist('Offset_from_UTC','var')
                            xlabel([datestr(mean(burst_times),'YYYY') ' (UTC+' Offset_from_UTC ')'])
                        else
                            xlabel(datestr(mean(burst_times),'YYYY'))
                        end
                        rotateXLabels( gca(), 90 )
                        h1.Color=[0 0 0];
                        %                 h2.Color=[ 0    0.4470    0.7410];
                        %                 h3.Color=[0.8500    0.3250    0.0980];
                        
                        fname=[whereput '/' instrument '_' location '_WAVE_PERIOD_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
                        % print('-dpng','-r200', fname);
                        % print('-dpdf', '-painters', fname);
                        %         export_fig(fname,'-pdf','-png','-transparent')
                        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
                        
                        
                end
                
                
            case 'zerocross'
                %              Hs: [100×1 double]
                %              Hz: [100×1 double]
                %              Tz: [100×1 double]
                %              Ts: [100×1 double]
                % Hmax
                
                f=figure;
                f.Position=[66 207 1600 450];
                f.Color='w';
                % set(gcf,'Position',[200 100 700 300])
                % set(gcf,'PaperPositionMode','auto')
                plot(burst_times,ocn.wave.Hs)
                hold on
                plot(burst_times,ocn.wave.Hz)
                plot(burst_times,ocn.wave.Hmax)
                legend('Hs','Hz','Hmax')
                ylabel('Wave height (m)')
                set(gca,'xlim',xlims,'xtick',xt)
                set(gca,'ylim',[0 ceil(max(ocn.wave.Hmax))])
                datetick('x','mm/dd','keeplimits','keepticks')
                rotateXLabels( gca(), 90 )
                xlabel(datestr(mean(burst_times),'YYYY'))
                title([instrument ' ' location ],'interpreter','none')
                fname=[whereput '/' instrument '_' location '_WAVE_HEIGHT_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
                % print('-dpng','-r200', fname);
                % print('-dpdf', '-painters', fname);
                %         export_fig(fname,'-pdf','-png','-transparent')
                savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
                
                f=figure;
                f.Position=[200 500 1600 450];
                f.Color='w';
                plot(burst_times,ocn.wave.Ts)
                hold on
                plot(burst_times,ocn.wave.Tz)
                set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
                legend('Ts','Tz')
                ylabel('Period (s)')
                title([instrument ' ' location ],'interpreter','none')
                datetick('x','mm/dd','keeplimits','keepticks')
                %         xlabel(datestr(mean(burst_times),'YYYY'))
                if exist('Offset_from_UTC','var')
                    xlabel([datestr(mean(burst_times),'YYYY') ' (UTC+' Offset_from_UTC ')'])
                else
                    xlabel(datestr(mean(burst_times),'YYYY'))
                end
                rotateXLabels( gca(), 90 )
                fname=[whereput '/' instrument '_' location '_WAVE_PERIOD_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
                % print('-dpng','-r200', fname);
                % print('-dpdf', '-painters', fname);
                %         export_fig(fname,'-pdf','-png','-transparent')
                savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
                
        end
        
        
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 subsample to less than 1 hr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Subsample to less than 1 hr...')


try
numdays=floor(ts(end)-ts(1)); % yh ?? -1
catch
numdays=floor(max(ts)-ts(1)); % yh ?? -1 yh test for stirling
end


clear inds spw tsub psub time_subsampled press_subsampled
sph=sample_rate_hz*3600;
co=0;
for sub_hr=1:numdays*24*(60/subsample_mins)
    co=co+1;
    spw=sph.*subsample_mins/60 ;  % samples per window
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
title([instrument ' ' location ' subsample interval = ' num2str(subsample_mins) ' minutes'],'interpreter','none');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create Hourly averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Create hourly averages...')

vec=datevec(ts(1));
v5  = vec(1,5)+vec(1,6)/60;
vec(5) = ceil(v5/30)*30;
vec(6) = 0;
Hstart=datenum(vec)+30/60/24;% this ensures correct number of samples at start

vec=datevec(ts(end));
v5  = vec(1,5)+vec(1,6)/60;
vec(5) = floor(v5/30)*30;
vec(6) = 0;
Hend=datenum(vec)-30/60/24; % this ensures correct number of samples at  end

Htime=[datenum(Hstart):1/24:datenum(Hend)];

Hpress=[]; clear Htemperature
for i=1:length(Htime);
    idx=find(ts>Htime(i)-30/60/24 & ts<Htime(i)+30/60/24);
    Hpress(i)=mean(p(idx));
    
    if exist('temp','var')
        Htemperature(i)=mean(temp(idx));
    end
end


% get avg pressure corresponding to burst averaging interval used in waves
switch calculate_waves
    case 'Y'
        burst_press=mean(ocn.wave.Burst_Data,2);
end



%% get IG waves (this overlaps with above)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch calculate_waves
    case 'Y'
        
%         [wavedata hs1  hm1 S]=get_chari_wave_analysis(ts,p,2); %[wavedata hs1  hm1 S]=get_chari_wave_analysis(time,data,srate); srate in Hz 2/second
        [wavedata hs1  hm1 S Pn Pstat]=get_chari_wave_analysis(ts,p,sample_rate_hz,burst_duration);
        
        
        
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
        title([instrument ' ' location],'interpreter','none');
        %         xlabel(datestr(mean(Htime),'YYYY'))
        if exist('Offset_from_UTC','var')
            xlabel([datestr(mean(ts),'YYYY') ' (UTC+' Offset_from_UTC ')'])
        else
            xlabel(datestr(mean(ts),'YYYY'))
        end
        % figure;
        %  plot(wavestata(1,:),wavestata(25,:))
        %  hold on
        %  xlabel('Time (days)')
        %  ylabel('IG height (m)')
        fname=[whereput '/' instrument '_' location '_IG_WAVES_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
        %         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
        
        
        
        % get hourly ig
        tig=wavedata(:,1);
        ig=wavedata(:,25);
        
        Hig=interp1(tig, ig,Htime);
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get lowpass sea level (lanczos or butterworth 38 hours)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    do t-tide to get tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(genpath('t_tide_v1'))
addpath(genpath(pwd))
% predict  tide & get residual from model data
clear names freq tidecon Htide

[names,freq,tidecon,Htide]=t_tide(Hpress,'interval',1,'start',Htime(1),'synthesis',1);
% mYout=t_predic(mts_sub,mnames,mfreq,mtidecon,'synthesis',1);
% calculate non-tidal residual
residual=Hpress-nanmean(Hpress)-Htide;

% get higher frequency tide output matching subsampled
% time_subsampled press_subsampled
tide_subsampled=t_predic(time_subsampled,names,freq,tidecon,'synthesis',1);
% calculate non-tidal residual
residual_subsampled=press_subsampled-nanmean(press_subsampled)-tide_subsampled;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % test
% tid=tide_subsampled+nanmean(press_subsampled);
% res=press_subsampled-tid;
% figure;
% hold on
% plot(time_subsampled,press_subsampled,'k')
% plot(time_subsampled,tid,'b')
% plot(time_subsampled,res,'r')
% legend('total','tide','residual','location','best')
% ylabel('Water level (m CD)')
% xlabel(['Date (' datestr(time_subsampled(1),'YYYY') ')'])
% set(gca,'xlim',xlims,'xtick',xt)
% datetick('x','mm/dd','keeplimits','keepticks')
% rotateXLabels( gca(), 90 )
% box on;
% grid on;
%
% f=figure;
% f.Position=[200 100 900 400];
% hold on
% plot(time_subsampled,press_subsampled,'k')
% plot(time_subsampled,tide_subsampled+nanmean(press_subsampled),'b')
% plot(time_subsampled,residual_subsampled,'r')
% legend('total','tide','residual','location','best')
% ylabel('Water level (m)')
% xlabel(['Date (' datestr(time_subsampled(1),'YYYY') ')'])
% set(gca,'xlim',xlims,'xtick',xt)
% datetick('x','mm/dd','keeplimits','keepticks')
% rotateXLabels( gca(), 90 )
% box on;
% grid on;
% title([instrument ' ' location ' Averaging interval = ' num2str(subsample_mins) ' minutes'],'interpreter','none');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          plot water levels
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch datum
    case 'CD'
        % plot chart datum
        f=figure;
        f.Position=[200 100 900 400];
        hold on
        plot(time_subsampled,press_subsampled,'k')
        plot(time_subsampled,tide_subsampled+mean(press_subsampled),'b')
        plot(time_subsampled,residual_subsampled,'r')
        % plot(time_subsampled(residual_subsampled<0),residual_subsampled(residual_subsampled<0),'g')
        legend('total','tide','residual','location','best')
        ylabel('Water level (m CD)')
        xlabel(['Date (' datestr(time_subsampled(1),'YYYY') ')'])
        yt=[floor(min(residual_subsampled)):.25:ceil(max(press_subsampled))];
        set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
        
        datetick('x','mm/dd','keeplimits','keepticks')
        rotateXLabels( gca(), 90 )
        box on;
        grid on;
        title([instrument ' ' location ' Averaging interval = ' num2str(subsample_mins) ' minutes'],'interpreter','none');
        
        fname=[whereput '/' instrument '_' location '_SEALEVEL_ChartDatum_' num2str(subsample_mins) 'mins_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
        
    case 'MSL'
        % plot msl
        f=figure;
        f.Position=[200 100 900 400];
        hold on
        plot(time_subsampled,press_subsampled-nanmean(press_subsampled),'k')
        plot(time_subsampled,tide_subsampled,'b')
        plot(time_subsampled,residual_subsampled,'r')
        legend('total','tide','residual','location','best')
        ylabel('Water level (m)')
        xlabel(['Date (' datestr(time_subsampled(1),'YYYY') ')'])
        yt=[floor(min(residual_subsampled)):.25:ceil(max(press_subsampled))];
        set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
        datetick('x','mm/dd','keeplimits','keepticks')
        rotateXLabels( gca(), 90 )
        box on;
        grid on;
        title([instrument ' ' location ' Averaging interval = ' num2str(subsample_mins) ' minutes'],'interpreter','none');
        
        fname=[whereput '/' instrument '_' location '_SEALEVEL_MSL_' num2str(subsample_mins) 'mins_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
        
        %                 h2.Color=[ 0    0.4470    0.7410];
        %                 h3.Color=[0.8500    0.3250    0.0980];
        
end

% plot hourly

switch datum
    case 'CD'
        % PLOT all sea levels
        f=figure;
        f.Color='w';
        set(gcf,'Position',[200 100 900 400])
        set(gcf,'PaperPositionMode','auto')
        hold on
        plot(Htime,Hpress,'k');
        plot(Htime,Htide+nanmean(Hpress),'b')
        plot(Htime,residual,'r')
        % plot(Htime,lp_wl-nanmean(lp_wl),'g')
        plot(Htime,lanczos_wl-nanmean(lanczos_wl),'g');
        % legend('Total  water level','Tide')
        ylabel('Sea level (m CD)')
        title([instrument ' ' location ],'interpreter','none')
        if exist('Offset_from_UTC','var')
            xlabel([datestr(mean(Htime),'YYYY') ' (UTC+' Offset_from_UTC ')'])
        else
            xlabel(datestr(mean(Htime),'YYYY'))
        end
        legend('Total','Tide','Residual','lowpass','location','best') % ,'lanczos lowpass'
        yt=[floor(min(residual)):.25:ceil(max(Hpress))];
        set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
        datetick('x','mm/dd','keeplimits','keepticks')
        rotateXLabels( gca(), 90 )
        box on;
        grid on;
        
        
        fname=[whereput '/' instrument '_' location '_SEALEVEL_ChartDatum_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
        %         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
        
end
%     case 'MSL'
% PLOT all sea levels
f=figure;
f.Color='w';
set(gcf,'Position',[200 100 900 400])
set(gcf,'PaperPositionMode','auto')
hold on
plot(Htime,Hpress-nanmean(Hpress),'k');
plot(Htime,Htide,'b')
plot(Htime,residual,'r')
% plot(Htime,lp_wl-nanmean(lp_wl),'g')
plot(Htime,lanczos_wl-nanmean(lanczos_wl),'g');
% legend('Total  water level','Tide')
ylabel('Sea level (m)')
title([instrument ' ' location ],'interpreter','none')
if exist('Offset_from_UTC','var')
    xlabel([datestr(mean(Htime),'YYYY') ' (UTC+' Offset_from_UTC ')'])
else
    xlabel(datestr(mean(Htime),'YYYY'))
end
legend('Total','Tide','Residual','lowpass','location','best') % ,'lanczos lowpass'
yt=[floor(min(residual)):.25:ceil(max(Hpress))];
set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
datetick('x','mm/dd','keeplimits','keepticks')
rotateXLabels( gca(), 90 )
box on;
grid on;


fname=[whereput '/' instrument '_' location '_SEALEVEL_MSL_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
% print('-dpng','-r200', fname);
% print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
% end

disp('T_tide analysis completed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot all [ need to edit for wave outputs]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            
            switch AnalysisMethod
                case 'spectral'
                    
                    switch SeparateSeaSwell
                        case 'yes'
                            
                            
                            ax1=subaxis(4,1,1)
                            plot(burst_times,ocn.wave.Hm0swell)
                            hold on
                            plot(burst_times,ocn.wave.Hm0sea)
                            plot(burst_times,ocn.wave.Hm0,'k')
                            legend('Hm0swell','Hm0sea','Hm0')
                            ylabel('Wave height (m)')
                            set(gca,'xlim',xlims,'xtick',xt)
                            set(gca,'ylim',[0 ceil(max(ocn.wave.Hm0))])
                            datetick('x','mm/dd','keeplimits','keepticks')
                            set(gca,'xticklabel',[])
                            %         rotateXLabels( gca(), 90 )
                            ax1.XTickLabelRotation=90;
                            %                             xlabel(datestr(mean(burst_times),'YYYY'))
                            title([instrument ' ' location ],'interpreter','none')
                            box on;
                            grid on;
                            
                            ax2=subaxis(4,1,2)
                            %                             h1=plot(burst_times,ocn.wave.Tp)
                            hold on
                            h2=plot(burst_times,ocn.wave.Tpswell)
                            h3=plot(burst_times,ocn.wave.Tpsea)
                            set(gca,'xlim',xlims,'xtick',xt,'ytick',[0:2:24])
                            %                             legend('Tp','swell','sea');
                            legend('swell','sea');
                            ylabel('Period (s)')
                            %                             title([instrument ' ' location ],'interpreter','none')
                            datetick('x','mm/dd','keeplimits','keepticks')
                            %         xlabel(datestr(mean(burst_times),'YYYY'))
                            set(gca,'xticklabel',[])
                            %                             if exist('Offset_from_UTC','var')
                            %                                 xlabel([datestr(mean(burst_times),'YYYY') ' (UTC+' Offset_from_UTC ')'])
                            %                             else
                            %                                 xlabel(datestr(mean(burst_times),'YYYY'))
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
                            plot(burst_times,ocn.wave.Hm0,'k')
                            legend('Hm0')
                            ylabel('Wave height (m)')
                            set(gca,'xlim',xlims,'xtick',xt)
                            set(gca,'ylim',[0 ceil(max(ocn.wave.Hm0))])
                            datetick('x','mm/dd','keeplimits','keepticks')
                            set(gca,'xticklabel',[])
                            %         rotateXLabels( gca(), 90 )
                            ax1.XTickLabelRotation=90;
                            %                             xlabel(datestr(mean(burst_times),'YYYY'))
                            title([instrument ' ' location ],'interpreter','none')
                            box on;
                            grid on;
                            
                            ax2=subaxis(4,1,2)
                            h1=plot(burst_times,ocn.wave.Tp)
                            set(gca,'xlim',xlims,'xtick',xt,'ytick',[0:2:24])
                            legend('Tp','swell','sea');
                            ylabel('Period (s)')
                            %                             title([instrument ' ' location ],'interpreter','none')
                            datetick('x','mm/dd','keeplimits','keepticks')
                            set(gca,'xticklabel',[])
                            %         xlabel(datestr(mean(burst_times),'YYYY'))
                            %                             if exist('Offset_from_UTC','var')
                            %                                 xlabel([datestr(mean(burst_times),'YYYY') ' (UTC+' Offset_from_UTC ')'])
                            %                             else
                            %                                 xlabel(datestr(mean(burst_times),'YYYY'))
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
                    plot(burst_times,ocn.wave.Hs)
                    hold on
                    plot(burst_times,ocn.wave.Hz)
                    plot(burst_times,ocn.wave.Hmax)
                    legend('Hs','Hz','Hmax')
                    ylabel('Wave height (m)')
                    set(gca,'xlim',xlims,'xtick',xt,'ytick',[0:0.5:5])
                    set(gca,'ylim',[0 ceil(max(ocn.wave.Hmax))])
                    datetick('x','mm/dd','keeplimits','keepticks')
                    %         rotateXLabels( gca(), 90 )
                    ax1.XTickLabelRotation=90;
                    set(gca,'xticklabel',[])
                    %                     xlabel(datestr(mean(burst_times),'YYYY'))
                    title([instrument ' ' location ],'interpreter','none')
                    box on;
                    grid on;
                    
                    ax2=subaxis(4,1,2)
                    plot(burst_times,ocn.wave.Ts)
                    hold on
                    plot(burst_times,ocn.wave.Tz)
                    set(gca,'xlim',xlims,'xtick',xt,'ytick',[0:2:24])
                    legend('Ts','Tz')
                    ylabel('Period (s)')
                    %                     title([instrument ' ' location ],'interpreter','none')
                    datetick('x','mm/dd','keeplimits','keepticks')
                    %         xlabel(datestr(mean(burst_times),'YYYY'))
                    set(gca,'xticklabel',[])
                    %                     if exist('Offset_from_UTC','var')
                    %                         xlabel([datestr(mean(burst_times),'YYYY') ' (UTC+' Offset_from_UTC ')'])
                    %                     else
                    %                         xlabel(datestr(mean(burst_times),'YYYY'))
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
                    datetick('x','mm/dd','keeplimits','keepticks')
                    set(gca,'xticklabel',[])
                    %         rotateXLabels( gca(), 90 )
                    ax3.XTickLabelRotation=90;
                    box on;
                    grid on;
                    
                case 'MSL'
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
                    yt=[floor(min(residual_subsampled)):.25:ceil(max(press_subsampled))];
                    set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
                    datetick('x','mm/dd','keeplimits','keepticks')
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
            plot(wavedata(:,1),wavedata(:,25))
            plot([Htime(1) Htime(end)],[0 0],'color',[.5 .5 .5],'linewidth',0.5)
            legend('IG waves','location','best')
            set(gca,'xlim',xlims,'xtick',xt)
            set(gca,'ytick',[-1:.05:1],'ylim',[0 .3])
            datetick('x','mm/dd','keeplimits','keepticks')
            %             rotateXLabels( gca(), 90 )
            ax3.XTickLabelRotation=90;
            ylabel('Height (m)')
            
            if exist('Offset_from_UTC','var')
                xlabel([datestr(mean(Htime),'YYYY') ' (UTC+' Offset_from_UTC ')'])
            else
                xlabel(datestr(mean(Htime),'YYYY'))
            end
            
            grid on
            
            
            fname=[whereput '/' instrument '_' location '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
            % print('-dpng','-r200', fname);
            % print('-dpdf', '-painters', fname);
            %         export_fig(fname,'-pdf','-png','-transparent')
            savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
            
        end
        
    case 'N'
        disp('NO waves plotted')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   interpolate waves to Htime and write to text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timestr=datestr(Htime','yyyy-mm-dd HH:MM:SS');
Time_WST=cellstr(timestr);



switch calculate_waves
    case 'Y'
        switch AnalysisMethod
            case 'spectral'
                
                switch SeparateSeaSwell
                    case 'yes'
                        
                        Hm0swell=interp1(burst_times,ocn.wave.Hm0swell,Htime);
                        Hm0sea=interp1(burst_times,ocn.wave.Hm0sea,Htime);
                        Hm0=interp1(burst_times,ocn.wave.Hm0,Htime);
                        Tp=interp1(burst_times,ocn.wave.Tp,Htime);
                        Tpswell=interp1(burst_times,ocn.wave.Tpswell,Htime);
                        Tpsea=interp1(burst_times,ocn.wave.Tpsea,Htime);
                        
                        Sealevel=reshape(Hpress,[],1);
                        Hm0swell=reshape(Hm0swell,[],1);
                        Hm0sea=reshape(Hm0sea,[],1);
                        Hm0=reshape(Hm0,[],1);
                        Hlowpass_wl=reshape(lowpass,[],1);
                        Htide=reshape(Htide,[],1);
                        Hresidual=reshape(residual,[],1);
                        Tp=reshape(Tp,[],1);
                        Tpswell=reshape(Tpswell,[],1);
                        Tpsea=reshape(Tpsea,[],1);
                        Hig=reshape(Hig,[],1);
                        
                        if exist('Htemperature','var')
                            T = table(timestr,Sealevel,Hlowpass_wl,Htide,Hresidual,Hm0,Tp,Hm0swell,Tpswell,Hm0sea,Tpsea,Hig,reshape(Htemperature,[],1));
                        else
                            T = table(timestr,Sealevel,Hlowpass_wl,Htide,Hresidual,Hm0,Tp,Hm0swell,Tpswell,Hm0sea,Tpsea,Hig);
                        end
                        
                        outf_tname=[whereput '/' instrument '_' location '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '_waves_hourly.csv'];
                        writetable(T,outf_tname,'Delimiter',',','QuoteStrings',true)
                        disp('Hourly wave and sea level data saved to csv text file')
                        
                    case 'no'
                        
                        Hm0=interp1(burst_times,ocn.wave.Hm0,Htime);
                        Tp =interp1(burst_times,ocn.wave.Tp,Htime);
                        
                        Sealevel=reshape(Hpress,[],1);
                        Hm0=reshape(Hm0,[],1);
                        Hlowpass_wl=reshape(lowpass,[],1);
                        Htide=reshape(Htide,[],1);
                        Hresidual=reshape(residual,[],1);
                        Tp=reshape(Tp,[],1);
                        Hig=reshape(Hig,[],1);
                        
                        if exist('Htemperature','var')
                            T = table(timestr,Sealevel,Hlowpass_wl,Htide,Hresidual,Hm0,Tp,Hig,reshape(Htemperature,[],1));
                        else
                            T = table(timestr,Sealevel,Hlowpass_wl,Htide,Hresidual,Hm0,Tp,Hig);
                        end
                        
                        outf_tname=[whereput '/' instrument '_' location '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '_waves_hourly.csv'];
                        writetable(T,outf_tname,'Delimiter',',','QuoteStrings',true)
                        disp('Hourly wave and sea level data saved to csv text file')
                end
                
            case 'zerocross'
                
                Hs =interp1(burst_times,ocn.wave.Hs,Htime);
                Hz =interp1(burst_times,ocn.wave.Hz,Htime);
                Hmax =interp1(burst_times,ocn.wave.Hmax,Htime);
                Ts =interp1(burst_times,ocn.wave.Ts,Htime);
                Tz =interp1(burst_times,ocn.wave.Tz,Htime);
                
                Sealevel=reshape(Hpress,[],1);
                Hs=reshape(Hs,[],1);
                Ts=reshape(Ts,[],1);
                Tz=reshape(Tz,[],1);
                Hz=reshape(Hz,[],1);
                Hmax=reshape(Hmax,[],1);
                Hlowpass_wl=reshape(lowpass,[],1);
                Htide=reshape(Htide,[],1);
                Hresidual=reshape(residual,[],1);
                Hig=reshape(Hig,[],1);
                
                if exist('Htemperature','var')
                    T = table(timestr,Sealevel,Hlowpass_wl,Htide,Hresidual,Hs,Ts,Hz,Tz,Hig,Htemperature');
                else
                    T = table(timestr,Sealevel,Hlowpass_wl,Htide,Hresidual,Hs,Ts,Hz,Tz,Hig);
                end
                
                outf_tname=[whereput '/' instrument '_' location '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '_waves_hourly.csv'];
                writetable(T,outf_tname,'Delimiter',',','QuoteStrings',true)
                
                disp('Hourly wave and sea level data saved to csv text file')
        end
        
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               save <1 hr sealevel data to text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sub_time=reshape(time_subsampled,[],1);
sub_tide=reshape(tide_subsampled,[],1);
sub_residual=reshape(residual_subsampled,[],1);
sub_press=reshape(press_subsampled,[],1);

if exist('temperature_subsampled','var')
%     sub_temp=temperature_subsampled;
    sub_temp=reshape(temperature_subsampled,[],1);    
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
outf_tname=[whereput '/' instrument '_' location '_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd') '_Sealevel_no_waves_' num2str(subsample_mins) '_minutes.csv'];
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
text(f12,ymiddle+260,' 12h','fontsize',[12])

% to add a vertical line at a period of 4 hours (f4)
f4=1/(4*60*60);
plot([f4 f4],ylim,'k:')
text(f4,ymiddle+100,' 4h','fontsize',[12])

% % 10 mins
f3=1/((15/60)*60*60);
plot([f3 f3],ylim,'k:')
text(f3,ytop,'15 min','fontsize',[12])

% f3=1/((1/60)*60*60);
% plot([f3 f3],ylim,'k:')
% text(f3,ytop,'1 min','fontsize',[12])

f20=1/300;
plot([f20 f20],ylim,'k:')
text(f20,ymiddle,'300 sec','fontsize',[12])
%
f20=1/30;
plot([f20 f20],ylim,'k:')
text(f20,ymiddle,'30 sec','fontsize',[12])
%
f12=1/15;
plot([f12 f12],ylim,'k:')
text(f12,.2,'15 sec','fontsize',[12])
title([instrument ' ' location],'interpreter','none')

f5=1/5;
plot([f5 f5],ylim,'k:')
text(f5,.6,'5 sec','fontsize',[12])
title([instrument ' ' location],'interpreter','none')

fname=[whereput '/' instrument '_' location '_SPECTRA_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
% print('-dpng','-r200', fname);
% print('-dpdf', '-painters', fname);
%         export_fig(fname,'-pdf','-png','-transparent')
savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use chari's function to calculate time-frequency plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
%         f24=log(1/(60*15));
%         plot(xlim,[f24 f24],'k-')
%         text(timef(100),f24+0.25,' 15 min','fontsize',[12])
        
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
        title([instrument ' ' location],'interpreter','none')
        set(gca,'xlim',[timef(1) timef(end)],'xtick',xt,'fontsize',font_size)
        datetick('x','mm/dd','keeplimits','keepticks')
        rotateXLabels( gca(), 90 )
        % set(gca,'xlim',[timef(1) timef(end)],'fontsize',font_size);
        
        
        
        fname=[whereput '/' instrument '_' location '_TIME_FREQ_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
        %         export_fig(fname,'-pdf','-png','-transparent')
        %         savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname) %this may cause matlab to hang up as figure too complicated
        switch use_export_fig
            case 'Y'
                export_fig(fname,'-png')
            case 'N'
                print(fname,'-dpng','-r300')
        end
    case 'N'
        disp('NO time-frequency plot created')
        
end

% 'VerticalAlignment', 'top'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% do highpass filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch do_highpass_filter
    case 'Y'
        disp('Applying highpass filter')
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
                    temp=i+datenum(str2num(datestr(mean(ts),'yyyy')),1,1);
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
                    title([instrument ' ' location ' highpass (<' num2str(max_period) ' s) sea level'],'interpreter','none')
                    
                    set(gca,'fontsize',12)
                    
                    
                    % pause
                    
                    fname=[whereput '/' instrument '_' location '_HIGHPASS_weekly_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
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
        title([instrument ' ' location ' highpass (<' num2str(max_period) ' s) sea level'],'interpreter','none')
        set(gca,'fontsize',12)
        box on
        grid on
        
        fname=[whereput '/' instrument '_' location '_HIGHPASS_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
        % print('-dpng','-r200', fname);
        % print('-dpdf', '-painters', fname);
        %         export_fig(fname,'-pdf','-png','-transparent')
        savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
        
    case 'N'
        disp('NO Highpass filter applied')
        
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the pressure corrected data to compare with raw [optional if created above]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % these are pressure corrected data
% bts=[burst_times(1):datenum(0,0,0,0,0,.5):burst_times(end)+datenum(0,0,0,0,59,59.5)]; %yh this is to read the pressure corrected raw data
% adj=reshape(ocn.wave.Eta',[1,281*7200]);
% 
% 
% % these are raw data highpass filtered
% 
%   figure;
%         hold on
%         set(gcf,'Position',[200 100 1200 400])
%         set(gcf,'PaperPositionMode','auto')
%         plot(bts,adj,'r');
%         plot(tsub,data,'b');
%         plot(wt(tf),wz(tf)-1,'r.')
%         set(gca,'xlim',xlims,'xtick',xt)
%         set(gca,'ylim',[floor(min(data)) ceil(max(data))],'ytick',[-1:.25:2])
%         datetick('x','mm/dd','keeplimits','keepticks')
%         rotateXLabels( gca(), 90 )
%         ylabel('Height (m)')
%         legend('pressure corrected','raw','location','best')
%         title([instrument ' ' location ' highpass (<' num2str(max_period) ' s) sea level'],'interpreter','none')
%         set(gca,'fontsize',12)
%         box on
%         grid on
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 save all data to mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch save_matfile
    case 'Y'
        % fname=[dir_ref instrument '_PROCESSED_DATA_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
        
        fname=[whereput '/' instrument '_' location '_PROCESSED_DATA_' datestr(ts(inds(1)) ,'yyyymmdd') '-' datestr(ts(inds(end)),'yyyymmdd')];
        save(fname)
        disp(['All processed variables saved to: ' fname '.mat'])
        
    case 'N'
        disp('NOT saving matfile')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

disp('DONE!!!!')
toc
