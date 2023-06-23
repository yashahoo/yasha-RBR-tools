% read list of sites to make urls & download recent bom obs for WA
% then concatinate and save to matlab structure
% yasha hetzel Jan 2017
%
% edits
% 2018-12-18 yh   line 152 temporarily excluded browse island to fix bug; need to
% recreate archive to include this site;
% 2019-01-23 yh fixed line 207 whwew typo caused script to fail resulting in
 % data gap from 19 dec 2018 to  23 jan 2019!
%  2019-3-19 added try catch to email error message
% 2019-10-15 modified to save each site .axf text file into its own directory because too many files
%            -ran move_data_to_directories.m to clean up old
%	     -fixed email password
% 2022-07-05 added witchcliff west in place of witchcliffe
clear all;
clc;


addpath(genpath('/home/yasha/matlab/m-files'));

%%

try

fname='WA_latest_BOM_obs_WEBLINKS_v2.csv'

[site_num,ident,Name,Lat,Lon,Height,link]=textread(fname,'%f %f %s %f %f %f %s','delimiter',',','headerlines',1);

% get link site id
link_id=nan(length(link),1);
for i=1:length(link)
temp=link{i};
link_id(i)=str2num(temp(end-8:end-4));
end

recent_dir='RECENT_DATA';
mkdir(recent_dir)


for i=1:length(link)
url=['http://www.bom.gov.au/fwo/IDW60801/IDW60801.' num2str(link_id(i)) '.axf']
% %url=['http://www.bom.gov.au/fwo/IDW60701/IDW60701.' num2str(link_id(i)) '.axf']
% fname=[recent_dir '/' sprintf('%02d',site_num(i)) '_wst_sta_' num2str(ident(i)) '_webID_' num2str(link_id(i)) '_' Name{i} '_' datestr(now,'yyyymmdd_HHMMSS')];
fname=[recent_dir '/' 'wst_sta_' num2str(ident(i)) '_webID_' num2str(link_id(i)) '_' Name{i} '_' datestr(now,'yyyymmdd_HHMMSS')]; % remove id at start of filename

try
websave(fname,url);
catch
disp(['INVALID URL: ' url])
end
end
disp('FInished downloading text files')

old_dir='OLD_DATA';
mkdir(old_dir)


%% read one of the data files (to be continued....)
recent_dir='RECENT_DATA';
% recent_dir='OLD_DATA';

D=dir(['./' recent_dir '/wst_sta*axf']);

% start loop here to  read all files
for si=1:length(D)
    
    
    
infile=[recent_dir '/' D(si).name]

% use function created with matlab data import 
% get first time to get only the lines with data
[wmo,name80,history_product80,local_date_time_full80,aifstime_utc80,lat,lon,apparent_t,delta_t,gust_kmh,gust_kt,air_temp,dewpt,press_msl,rain_trace80,rel_hum,wind_dir80,wind_spd_kmh,wind_spd_kt] = importBOM_latest(infile);

dud=~isnan(wmo);
goodi=find(dud>0);

if ~isempty(goodi)

startline=goodi(1);
endline=goodi(end);


keep infile startline endline si recent_dir old_dir D recent site_num Name ident

% get the good data
[wmo,name80,history_product80,local_date_time_full80,aifstime_utc80,lat,lon,apparent_t,delta_t,gust_kmh,gust_kt,air_temp,dewpt,press_msl,rain_trace80,rel_hum,wind_dir80,wind_spd_kmh,wind_spd_kt] = importBOM_latest(infile,startline,endline);

% convert time into matlab serial time
local_date_time_full80=strrep( local_date_time_full80(:,1),'"',''); % get rid of starnge double quotes
aifstime_utc80=strrep( aifstime_utc80(:,1),'"','');
mtime_UTC=datenum(aifstime_utc80,'yyyymmddHHMMSS');
mtime_local=datenum(local_date_time_full80,'yyyymmddHHMMSS');


% sort the data ascending
[dud idx]=sort(mtime_UTC);

mtime_UTC=mtime_UTC(idx);
mtime_local=mtime_local(idx);
wmo = wmo(idx);
name80 = name80(idx); name80 = strrep( name80(:),'"','');
history_product80 =  history_product80(idx); history_product80 = strrep(history_product80(:),'"','');
local_date_time_full80 = local_date_time_full80(idx); 
aifstime_utc80= aifstime_utc80(idx);
lat= lat(idx);
lon= lon(idx);
apparent_t= apparent_t(idx);
delta_t= delta_t(idx);
gust_kmh= gust_kmh(idx);
gust_kt= gust_kt(idx);
air_temp= air_temp(idx);
dewpt = dewpt(idx);
press_msl =press_msl(idx);
rain_trace80 = rain_trace80(idx);
rel_hum = rel_hum(idx);
wind_dir80 = wind_dir80(idx);wind_dir80 = strrep( wind_dir80(:),'"','');
wind_spd_kmh = wind_spd_kmh(idx);
wind_spd_kt = wind_spd_kt(idx);

disp(['reading ' name80{1}])

% write recent data to structure
recent(si).name=name80;
recent(si).BOM_station_ID=history_product80;
recent(si).wmo_id=wmo;
recent(si).mtime_UTC=mtime_UTC;
recent(si).lat=lat;
recent(si).lon=lon;
recent(si).apparent_t=apparent_t;
recent(si).delta_t=delta_t;
recent(si).gust_kmh=gust_kmh;
recent(si).gust_kt=gust_kt;
recent(si).air_temp=air_temp;
recent(si).dewpt=dewpt;
recent(si).press_msl=press_msl;
recent(si).rain_trace=rain_trace80;
recent(si).rel_hum=rel_hum;
recent(si).wind_dir=wind_dir80;
recent(si).wind_spd_kmh=wind_spd_kmh;
recent(si).wind_spd_kt=wind_spd_kt;
% recent(si).sitelist=Name;
% recent(si).site_num=site_num;

else
 disp(['NO GOOD DATA IN: ' infile])
end
end

 %%
% move recent data to old folders
for i=1:length(Name)
disp(['MOVING DATA.... i= ' num2str(i)])
try
% %sitefilename=['./RECENT_DATA/*.axf'
%  %movefile('./RECENT_DATA/*.axf','./OLD_DATA')
% sitefilename=['./RECENT_DATA/' sprintf('%02d',site_num(i)) '_wst_sta_' num2str(ident(i)) '*_' Name{i} '*.axf'];
% sitedir=[old_dir '/' sprintf('%02d',site_num(i)) '_wst_sta_' num2str(ident(i)) '_' Name{i}];

sitefilename=['./RECENT_DATA/wst_sta_' num2str(ident(i)) '*_' Name{i} '*.axf'];
sitedir=[old_dir '/' 'wst_sta_' num2str(ident(i)) '_' Name{i}];

%mkdir(sitedir);
disp([sitefilename ' --> ' sitedir '/'])
 movefile(sitefilename,[sitedir '/'])

 %disp('Recent data moved to OLD')
  
catch
disp(['ERROR...........' sitefilename ' --> ' sitedir '/' '.....OR....may not exist'])
   % disp('No new data to MOVE!!')
end
end




%archive=recent; 
%save BOM_obs_archived archive  Name site_num % only for first time
%% Combine older and more recent data from structure
load BOM_obs_archived
keep recent archive Name site_num

% get list of names in archive
archnames=[];
for ai=1:length(archive)
    try
        archnames{ai}=archive(ai).name{1};
    catch
        archnames{ai}='NaN';
        archive(ai).name{1}={'NaN'};
    end
end

% get list of names in recent
recnames=[];
for ri=1:length(recent)
    recnames{ri}=recent(ri).name{1};
end

% find which sites are different
[C1,IA1,in_archive_not_recent] = union(recnames,archnames);
[C2,IA2,in_recent_not_archive] = union(archnames,recnames);
for i=1:length(in_archive_not_recent)
disp([archive(in_archive_not_recent(i)).name{1} ' exists in archive(' num2str(in_archive_not_recent(i)) ') but not recent'])
end
for i=1:length(in_recent_not_archive)
disp([recent(in_recent_not_archive(i)).name{1} ' exists in recent(' num2str(in_recent_not_archive(i)) ') but not archive'])
end

disp('Combining OLD and NEW data in structure')
co_miss_data=0;
for si=1:length(recent)%-1 % loop through sites -1 added to temporarily fix prblime where bvrowse island data are missing for earlier archive
    
    %     find if recent site exists in archive and get index for it
    tmpname=recent(si).name{1};
    
    asi=find(strcmpi(tmpname,archnames)) ;
    if isempty(asi)
        co_miss_data=co_miss_data+1; % counter of missing data to not overwrite old data
        asi=length(archive)+co_miss_data;
        old_data_exists='N';
    else
        old_data_exists='Y';
    end
    
    disp([num2str(si) '. ' tmpname ' old data exists: '  old_data_exists])
    
    switch old_data_exists
        case 'Y'
            % ---------------------------------------------------------------------- %
            % ------------------------- if old data exist ------------------------- %
            % ---------------------------------------------------------------------- %
            try
                atime=archive(asi).mtime_UTC;
                if isempty(atime)
                    disp(['atime is empty for site: ' num2str(si) ' ' tmpname ])
                end
                ltime=atime(end);
                
                
                rti =find(recent(si).mtime_UTC>ltime);
                disp(['At site: ' recent(si).name{1} '  new data added: '])
%                 disp(datestr(recent(si).mtime_UTC(rti)))
                disp([datestr(recent(si).mtime_UTC(rti(1))) ' to ' datestr(recent(si).mtime_UTC(rti(end)))])
                
                %try
                dtgap=(recent(si).mtime_UTC(rti(1))-ltime)*24;
                if dtgap>0.51
                    disp(['Data gap of ' num2str(dtgap) ' hrs or ' num2str(dtgap./24) ' days'])
                end
                
                
                
                f = fieldnames(recent);
                for i = 1:length(f)
                    combined(si).(f{i}) = [archive(asi).(f{i});recent(si).(f{i})(rti)];
                end
                
                
                
                
            catch
                disp(['NO NEW DATA AVAILABLE FOR SITE: '  num2str(si) ' ' tmpname   ])
                f = fieldnames(recent);
                for i = 1:length(f)
                    combined(si).(f{i}) = archive(asi).(f{i});
                end
            end
            % ---------------------------------------------------------------------- %
            % ----------------------- if NO old data exist ------------------------- %
            % ---------------------------------------------------------------------- %
        case 'N'
            
            try
                % atime=archive(asi).mtime_UTC;
                atime=[];
                if isempty(atime)
                    disp(['atime is empty for site: ' num2str(si) ' ' tmpname ])
                end
                % ltime=atime(end);
                ltime=recent(si).mtime_UTC(1)-1/3600/24; % get 1 sec before first recent time
                
                
                rti =find(recent(si).mtime_UTC>ltime);
                disp(['At site: ' recent(si).name{1} '  new data added: '])
%                 disp(datestr(recent(si).mtime_UTC(rti)))
                disp([datestr(recent(si).mtime_UTC(rti(1))) ' to ' datestr(recent(si).mtime_UTC(rti(end)))])
                
                %try
                dtgap=(recent(si).mtime_UTC(rti(1))-ltime)*24;
                if dtgap>0.51
                    disp(['Data gap of ' num2str(dtgap) ' hrs or ' num2str(dtgap./24) ' days'])
                else
                    disp(['No archive data available... adding only new...'])
                end
                
                
                
                f = fieldnames(recent);
                for i = 1:length(f)
                    %     combined(si).(f{i}) = [archive(asi).(f{i});recent(si).(f{i})(rti)];
                    combined(asi).(f{i}) = [recent(si).(f{i})(rti)];
                    
                end
                
                
                
                
            catch
                disp(['NO NEW DATA AVAILABLE FOR SITE: '  num2str(si) ' ' tmpname   ])
                f = fieldnames(recent);
                for i = 1:length(f)
                    combined(si).(f{i}) = archive(asi).(f{i});
                end
            end
            
    end
    % ---------------------------------------------------------------------- %
    
end


% throw out any empty stations
% find indexes in structure
co=0;empty_siti=[];
for i=1:length(combined);
if isempty(combined(i).name)
    co=co+1;
empty_siti(co)=i;
end
end
% delete them
if ~isempty(empty_siti)
combined(empty_siti)=[];
end

% make new site list for combined
% site_list=Name;
co=0;site_list=[];site_num=[];
for i=1:length(combined);
    co=co+1;
    site_list{co}=combined(i).name{1};
    site_num(co)=i;
end
    


try
    archive=combined;site_num=site_num';
    save BOM_obs_archived archive site_num site_list
    disp('Data updated and saved in archive...now calculating u,v etc. and saving again (slow)')
catch
    disp('No new data to combine!!')
end

% %% copy data to public_html
% copyfile('BOM_obs_archived.mat', '/home/yasha/public_html/')
% copyfile('SA_BOM_obs_archived.mat', '/home/yasha/public_html/')

 %% test plot
 close all
 
 %for si=1:14% 32
 
 %t=archive(si).mtime_UTC;
 %dat=archive(si).wind_spd_kt;
 %dat(dat<-1000)=nan;
% dat2=archive(si).gust_kt;
 %dat2(dat2<-1000)=nan;
 
 %figure;
 %plot(t+8/24,dat,'k')
 %hold on
 %plot(t+8/24,dat2,'b')
 %datetick('x','dd-mmm HH:MM','keeplimits')
 %ylabel('Wind speed (kts)')
 %xlabel('Date')
 %title([archive(si).name{1} '  most recent data: ' datestr(t(end)+8/24) ' wst'])
 
 %end

 
 %%
 
%  get direction in degrees from  wind direction text
 
%  get u,v
dirs=[0:22.5:337.5 0 0];
missingdir=[0:22.5:337.5 0 NaN];
dirstr={'N' 'NNE' 'NE' 'ENE' 'E' 'ESE' 'SE' 'SSE' 'S' 'SSW' 'SW' 'WSW' 'W' 'WNW' 'NW' 'NNW' 'CALM'};

for si=1:length(archive)
    for i=1:length(archive(si).mtime_UTC)

try
archive(si).wind_dir_deg(i)=dirs(strcmp(archive(si).wind_dir(i),dirstr));
catch
archive(si).wind_dir_deg(i)=NaN;
end

         %make missing data be nans
        archive(si).wind_spd_kmh(archive(si).wind_spd_kmh<0)=NaN;

        archive(si).wind_spd_ms(i)=archive(si).wind_spd_kmh(i).*1000./3600;
        [u v]=(compass2cart(archive(si).wind_dir_deg(i),archive(si).wind_spd_ms(i)));
        try
        archive(si).u(i)=-u;
        archive(si).v(i)=-v;
        catch
           archive(si).u(i)= NaN;
           archive(si).v(i)= NaN;
        end
        
         %make missing data be nans
        archive(si).gust_kmh(archive(si).gust_kmh<0)=NaN;
        
        archive(si).gust_ms(i)=archive(si).gust_kmh(i).*1000./3600;
        [ug vg]=(compass2cart(archive(si).wind_dir_deg(i),archive(si).gust_ms(i)));
        try
        archive(si).ug(i)=-ug;
        archive(si).vg(i)=-vg;
        catch
           archive(si).u(i)= NaN;
           archive(si).v(i)= NaN;
        end

%         % this is to make accurate directions but it looks bad to have gaps
%         if isnan(missingdir(strcmp(archive(si).wind_dir(i),dirstr)))
%             archive(si).u(i)=NaN;
%             archive(si).v(i)=NaN;
%             archive(si).ug(i)=NaN;
%             archive(si).vg(i)=NaN;
%             disp(['making dir NaN for site: '  num2str(si)   '  ' datestr(archive(si).mtime_UTC(i) )])
%         end
        
        
        
        
        clear u v
    end
end

% 
% save BOM_obs_archived archive site_num site_list
% 
% disp('Data updated and saved in archive with U V components and dir degrees')

%% plot wind arrow plot
% for si=7:10
%     
% %WindArrows3(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,datestr(floor(archive(si).mtime_UTC(end))-4),5,1,15)
% WindArrows4(archive(si).u,archive(si).v,archive(si).mtime_UTC+8/24,datestr(floor(archive(si).mtime_UTC(end))-4),6,1,12,15)
% %title([archive(si).name{1} ])
% title([archive(si).name{1} '  most recent data: ' datestr(archive(si).mtime_UTC(end)+8/24) 'wst'])
% end

%%close all;

whenrunlast=now;

save BOM_obs_archived archive site_num site_list whenrunlast
disp(['Data updated and saved in archive including u and v components at: ' datestr(whenrunlast)]) 


catch ME
  
[fullFileName]=WarnUser(ME)
   
% email broken!!
% %     email if there is a problem
% % do this once
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','E_mail','yasha.hetzel@gmail.com');
% setpref('Internet','SMTP_Username','yasha.hetzel');
% setpref('Internet','SMTP_Password','Tombies9253');
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');
% 
% 
% 
% % then use this
% sendmail('yashahoo@gmail.com','Error downloading BOM files',['Error downloading Bom data on ' datestr(now) ' '  ME.message],fullFileName)
% %keyboard


end

%% copy data to public_html
%copyfile('BOM_obs_archived.mat', '/home/yasha/public_html/')
copyfile('BOM_obs_archived.mat','/home2/yasha/BOM/')
%copyfile('SA_BOM_obs_archived.mat', '/home/yasha/public_html/')
disp('FINISHED COPYING to yasha@rottnest:/home2/yasha/BOM/')
