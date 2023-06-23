function [archive]=SA_plot_date_BOM_wind(sites,date,numdays,ymax)

% plot wind arrow plots for recent bom obs for SA (2017-03-15 to present)
% **need to rebuild data file to fix missing data, but ok for now.
% data stored in matlab structurein BOM_obs_archived.mat
% USAGE:
% [archive]=SA_plot_date_BOM_wind(sites,date,numdays,ymax);
% site =  station number (can be vector)
% date in string format (starting)
% numdays= number of days to plot (best under 10)
% ymax= ymax in m/s
% example: [archive]=SA_plot_date_BOM_wind([7:10],'2018-01-27',6,20) 
%  [13:15] 
% sites is site number vector or single number corresponding to 
% first column in SA_latest_BOM_obs_WEBLINKS.csv  
% 
% 1     18106	NULLARBOR	-31.45	130.9	64
% 2     18200	THEVENARD	-32.18	133.63	8.54
% 3     18191	COLES_POINT	-34.37	135.37	28
% 4     18192	PORT_LINCOLN AWS	-34.6	135.88	8.5
% 5     18115	NEPTUNE_ISLAND	-35.34	136.12	32
% 6     18120	WHYALLA_AERO	-33.05	137.52	9.3
% 7     18201	PORT_AUGUSTA_AERO	-32.51	137.72	13.97
% 8     21139	PORT_PIRIE_AERODROME AWS	-33.24	138	12
% 9     22049	STENHOUSE_BAY	-35.28	136.94	42
% 10	22823	CAPE_BORDA	-35.75	136.6	158
% 11	22803	CAPE_WILLOUGHBY	-35.84	138.13	55
% 12	23052	BLACK_POLE	-34.73	138.47	6.27
% 13	23034	ADELAIDE_AIRPORT	-34.95	138.52	2
% 14	26095	THE_LIMESTONE	-36.97	139.72	17
% 
% yasha hetzel Jan 2017
% addpath(genpath('/home/yasha/matlab/m-files'));
keep sites date numdays ymax;
% close all
% clc;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLot wind arrows for last  5 days


% % load data
% load SA_BOM_obs_archived.mat

% load data
if exist('SA_BOM_obs_archived.mat','file')
    if ismember('archive',evalin('base','who'))
        archive=evalin('base','archive');
    elseif ~exist('archive','var')       
        disp('Loading available data file...')
        load SA_BOM_obs_archived.mat
    end
    if datenum(date)>archive(sites(1)).mtime_UTC(end)+8/24
        disp('Current data file not available...downloading (must be connected to UWA network)')
        [datafile]=websave('SA_BOM_obs_archived.mat','http://130.95.29.59/~yasha/SA_BOM_obs_archived.mat');
        disp(['Download complete, most recent data available: ' datestr(archive(sites(1)).mtime_UTC(end)+8/24) 'WST'])
        disp('Loading data file...')
        load(datafile)
        
    end
    
elseif ~exist('SA_BOM_obs_archived.mat','file')
    
    try
        disp('Data file not available...downloading (must be connected to UWA network)')
        [datafile]=websave('SA_BOM_obs_archived.mat','http://130.95.29.59/~yasha/SA_BOM_obs_archived.mat');
        disp(['Download complete, most recent data available: ' datestr(archive(sites(1)).mtime_UTC(end)+8/24) ' WST'])
        disp('Loading data file...')
        load(datafile)
        
    catch
        f = errordlg('Unable to download file...check you are on UWA network','Error');
        return;
    end
    
end

if datenum(date)>archive(sites(1)).mtime_UTC(end)+8/24 || datenum(date)<archive(sites(1)).mtime_UTC(1)+8/24
    %             disp('Requested Date range not available...')
    f = errordlg('Requested Date range not available! Must be between 2017-01-18 and 6am today','Error');
    return;
end


%%
% plot wind arrow plot
for si=sites
    disp(['Plotting wind arrows for site: ' num2str(si)])
%%WindArrows3(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,datestr(floor(archive(si).mtime_UTC(end))-4),5,1,15)
WindArrows4(archive(si).u,archive(si).v,archive(si).mtime_UTC+9.5/24,date,numdays,1,12,ymax)

%WindArrows4(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,datestr(floor(archive(si).mtime_UTC(end))-6),8,1,12,15)

%%title([archive(si).name{1} ])
 title([archive(si).name{1} ])
xlabel('Date (ACT)')
set(gcf,'color','w');
end

end




