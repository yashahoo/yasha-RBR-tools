function [archive]=plot_date_BOM_wind_v2(sites,date,numdays,ymax,ymax2)
% Plot wind arrow plots for recent bom obs for WA (2017-01-17 to present)
% data stored in matlab structurein BOM_obs_archived.mat
% USAGE:
% [archive]=plot_date_BOM_wind_v2(sites,date,numdays,ymax,ymax2)
% site =  station number (can be vector)
% date in string format (starting) in WST (UTC+8)
% numdays= number of days to plot (best under 10)
% ymax= ymax in m/s
% example: [archive]=plot_date_BOM_wind_v2([7:10],'2018-01-27',6,20,20)  plots sites 7 to 10 (Cape Leeuwin to Busselton)
%  [13:15] 
% for recent data plots: datestr(floor(now)-6)
% sites is site number vector or single number corresponding to 
% first column in WA_latest_BOM_obs_WEBLINKS.csv  (for SW use [7:10]
% This functions call the function WindArrows4.m which can be used to plot arrows for any time,u,v data
% e.g. WindArrows4(u,v,wtime,begintime,numdays,interval,xtickhrs,ymax)
% WindArrows6(u,v,wtime,begintime,numdays,interval,xtickhrs,ymax,ymax2)
% http://www.bom.gov.au/fwo/IDW60801/IDW60801.94600.axf
% yasha hetzel Jan 2017
% 
% sites:
% site_num	ID	Name	Lat	Lon	Height
% 1. IDW60801    Witchcliffe West       -34 S 115.1 E
% 2. IDW60801    Eucla       -31.7 S 128.9 E
% 3. IDW60801    Red Rocks Point       -32.2 S 127.5 E
% 4. IDW60801    Rowley Shoals       -17.5 S 119 E
% 5. IDW60801    Adele Island       -15.5 S 123.2 E
% 6. IDW60801    Browse Island       -14.1 S 123.5 E
% 7. IDW60801    Broome       -17.9 S 122.2 E
% 8. IDW60801    Port Hedland       -20.4 S 118.6 E
% 9. IDW60801    Karratha       -20.7 S 116.8 E
% 10. IDW60801    Bedout Island       -19.6 S 119.1 E
% 11. IDW60801    Learmonth       -22.2 S 114.1 E
% 12. IDW60801    Onslow Airport       -21.7 S 115.1 E
% 13. IDW60801    Barrow Island       -20.9 S 115.4 E
% 14. IDW60801    Carnarvon       -24.9 S 113.7 E
% 15. IDW60801    Shark Bay Airport       -25.9 S 113.6 E
% 16. IDW60801    North Island       -28.3 S 113.6 E
% 17. IDW60801    Geraldton Airport       -28.8 S 114.7 E
% 18. IDW60801    Perth Airport       -31.9 S 116 E
% 19. IDW60801    Rottnest Island       -32 S 115.5 E
% 20. IDW60801    Ocean Reef       -31.8 S 115.7 E
% 21. IDW60801    Swanbourne       -32 S 115.8 E
% 22. IDW60801    Garden Island       -32.2 S 115.7 E
% 23. IDW60801    Hillarys Point Boat Harbour       -31.8 S 115.7 E
% 24. IDW60801    Cape Leeuwin       -34.4 S 115.1 E
% 25. IDW60801    Cape Naturaliste       -33.5 S 115 E
% 26. IDW60801    Busselton Airport       -33.7 S 115.4 E
% 27. IDW60801    Esperance       -33.8 S 121.9 E
% 28. IDW60801    Hopetoun North       -33.9 S 120.1 E
% 29. IDW60801    Bunbury       -33.4 S 115.6 E
% 30. IDW60801    Mandurah       -32.5 S 115.7 E
% 31. IDW60801    Walpole North       -34.9 S 116.7 E
% 32. IDW60801    Albany Airport       -34.9 S 117.8 E
% 
%  updates:
% 20190320 - added check for file availablity and download if not up to
% date enough; also to not reload if archive variable already loaded into worksspace
% 20200128 - use new v6 wind plotting function to adjust axis limits
% 20220706 - 

%addpath(genpath('/home/yasha/matlab/m-files'));
keep sites date numdays ymax ymax2;
% close all
% clc;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLot wind arrows 


% load data
if exist('BOM_obs_archived.mat','file')
    if ismember('archive',evalin('base','who'))
        archive=evalin('base','archive');
    elseif ~exist('archive','var')       
        disp('Loading available data file...')
        load BOM_obs_archived.mat
    end
    if datenum(date)>archive(sites(1)).mtime_UTC(end)+8/24
        disp('Current data file not available...downloading from https://sealevelx.ems.uwa.edu.au/plot/yasha/')
%         [datafile]=websave('BOM_obs_archived.mat','http://130.95.29.59/~yasha/BOM_obs_archived.mat');
        [datafile]=websave('BOM_obs_archived.mat','https://sealevelx.ems.uwa.edu.au/plot/yasha/BOM_obs_archived.mat');
        disp(['Download complete, most recent data available: ' datestr(archive(sites(1)).mtime_UTC(end)+8/24) 'WST'])
        disp('Loading data file...')
        load(datafile)
        
    end
    
elseif ~exist('BOM_obs_archived.mat','file')
    
    try
        disp('Data file not available...downloading from https://sealevelx.ems.uwa.edu.au/plot/yasha/')
%         [datafile]=websave('BOM_obs_archived.mat','http://130.95.29.59/~yasha/BOM_obs_archived.mat');
        [datafile]=websave('BOM_obs_archived.mat','https://sealevelx.ems.uwa.edu.au/plot/yasha/BOM_obs_archived.mat');
        
        disp('Loading data file...')
        load(datafile)
        disp(['Download complete, most recent data available: ' datestr(archive(sites(1)).mtime_UTC(end)+8/24) ' WST'])
        
    catch
        f = errordlg('Unable to download or load file...','Error');
        return;
    end
    
end

if datenum(date)>archive(sites(1)).mtime_UTC(end)+8/24 || datenum(date)<archive(sites(1)).mtime_UTC(1)+8/24
    %             disp('Requested Date range not available...')
%     f = errordlg('Requested Date range not available! Must be between 2017-01-18 and 6am today','Error');
     f = errordlg(['Requested Date range not available! Must be between ' datestr(archive(sites(1)).mtime_UTC(1)+8/24,'yyyy-mm-dd HH:MM') ' and ' datestr(archive(sites(1)).mtime_UTC(end)+8/24,'yyyy-mm-dd HH:MM')],'Error');
    return;
end

%%
% plot wind arrow plot
for si=sites
    disp(['Plotting wind arrows for site: ' num2str(si)])
% %%WindArrows3(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,datestr(floor(archive(si).mtime_UTC(end))-4),5,1,15)
% WindArrows4(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,date,numdays,1,12,ymax)
figure;
% WindArrows6(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,date,numdays,1,24,ymax,ymax2)
WindArrows6(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,date,numdays,1,24,ymax,ymax2)
% %WindArrows4(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,datestr(floor(archive(si).mtime_UTC(end))-6),8,1,12,15)

%%title([archive(si).name{1} ])
 title([archive(si).name{1} ])
xlabel(['Date ' datestr(datenum(date),'yyyy') ' (AWST)'])
set(gcf,'color','w');
f=gcf;
f.Position=[127 734 2175 484];
ylabel('Windspeed (m s^-^1)')
end


end




