function [archive]=plot_date_BOMwind(sites,date,numdays,ymax)
% Plot wind arrow plots for recent bom obs for WA (2017-01-17 to present)
% data stored in matlab structurein BOM_obs_archived.mat
% USAGE:
% [archive]=plot_date_BOMwind(sites,date,numdays,ymax)
% site =  station number (can be vector)
% date in string format (starting) in WST (UTC+8)
% numdays= number of days to plot (best under 10)
% ymax= ymax in m/s
% example: [archive]=plot_date_BOM_wind([7:10],'2018-01-27',6,20)  plots sites 7 to 10 (Cape Leeuwin to Busselton)
%  [13:15] 
% sites is site number vector or single number corresponding to 
% first column in WA_latest_BOM_obs_WEBLINKS.csv  (for SW use [7:10]
% This functions call the function WindArrows4.m which can be used to plot arrows for any time,u,v data
% e.g. WindArrows4(u,v,wtime,begintime,numdays,interval,xtickhrs,ymax)
% http://www.bom.gov.au/fwo/IDW60801/IDW60801.94600.axf
% yasha hetzel Jan 2017
%
% sites:
% site_num	ID	Name	Lat	Lon	Height
% 1 	11003	EUCLA           -31.68	128.9	93.1
% 2     11053	RED_ROCKS_POINT	-32.2	127.53	3.1
% 3     9789	ESPERANCE       -33.83	121.89	25
% 3     9961	HOPETOUN_NORTH	-33.93	120.13	26
% 5     9999	ALBANY_AIRPORT	-34.94	117.82	68.4
% 6     9998	NORTH_WALPOLE	-34.95	116.72	73
% 7     9518	CAPE_LEEUWIN	-34.37	115.14	13
% 8     9746	WITCHCLIFFE     -34.03	115.1	80
% 9     9519	CAPE_NATURALISTE-33.54	115.02	109
% 10	9603	BUSSELTON_AERO	-33.69	115.4	16.3
% 11	9965	BUNBURY         -33.36	115.64	5.2
% 12	9977	MANDURAH        -32.52	115.71	3
% 13	9256	GARDEN_ISLAND_HSF-32.24	115.68	6
% 14	9193	ROTTNEST_ISLAND	-32.01	115.5	43.1
% 15	9021	PERTH_AIRPORT	-31.93	115.98	15.4
% 16	9215	SWANBOURNE      -31.96	115.76	40.96
% 17	9265	HILLARYS_BOAT_HARBOUR_NTC_AWS	-31.83	115.74	0
% 18	9214	OCEAN_REEF      -31.76	115.73	10
% 19	8315	GERALDTON_AIRPORT	-28.8	114.7	29.7
% 20	8290	NORTH_ISLAND	-28.3	113.59	2.2
% 21	6105	SHARK_BAY_AIRPORT	-25.89	113.58	33.8
% 22	6011	CARNARVON_AIRPORT	-24.89	113.67	4
% 23	5007	LEARMONTH_AIRPORT	-22.24	114.1	5
% 24	5017	ONSLOW_AIRPORT	-21.67	115.11	10.5
% 25	5094	BARROW_ISLAND_AIRPORT	-20.87	115.41	6.4
% 26	4083	KARRATHA_AERO	-20.71	116.77	5.3
% 27	4032	PORT_HEDLAND_AIRPORT	-20.37	118.63	6.4
% 28	4100	BEDOUT_ISLAND	-19.59	119.1	8.6
% 29	3003	BROOME_AIRPORT	-17.95	122.24	7.42
% 30	200713	ROWLEY_SHOALS	-17.52	118.95	0
% 31	200735	ADELE_ISLAND	-15.51	123.16	4.6
% 32	200784	BROWSE_ISLAND	-14.11	123.55	3.7
% 
%  updates:
% 20190320 - added check for file availablity and download if not up to
% date enough; also to not reload if archive variable already loaded into
% worksspace

%addpath(genpath('/home/yasha/matlab/m-files'));
keep sites date numdays ymax;
% close all
% clc;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLot wind arrows for last  5 days


% load data
if exist('BOM_obs_archived.mat','file')
    if ismember('archive',evalin('base','who'))
        archive=evalin('base','archive');
    elseif ~exist('archive','var')       
        disp('Loading available data file...')
        load BOM_obs_archived.mat
    end
    if datenum(date)>archive(sites(1)).mtime_UTC(end)+8/24
        disp('Current data file not available...downloading (must be connected to UWA network)')
        [datafile]=websave('BOM_obs_archived.mat','http://130.95.29.59/~yasha/BOM_obs_archived.mat');
        disp(['Download complete, most recent data available: ' datestr(archive(sites(1)).mtime_UTC(end)+8/24) 'WST'])
        disp('Loading data file...')
        load(datafile)
        
    end
    
elseif ~exist('BOM_obs_archived.mat','file')
    
    try
        disp('Data file not available...downloading (must be connected to UWA network)')
        [datafile]=websave('BOM_obs_archived.mat','http://130.95.29.59/~yasha/BOM_obs_archived.mat');
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
WindArrows4(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,date,numdays,1,12,ymax)

%WindArrows4(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,datestr(floor(archive(si).mtime_UTC(end))-6),8,1,12,15)

%%title([archive(si).name{1} ])
 title([archive(si).name{1} ])
xlabel('Date (AWST)')
set(gcf,'color','w');
end


end




