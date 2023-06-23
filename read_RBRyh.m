function [netcdf_out]=read_RBRyh(filename,location,whereput) %,mslp_time,mslp)
% this reads RBR files using ruskin toolbox rbr-rsktools
% it saves a netcdf and textfile of the data
% usage:
% [netcdf_out ]=read_RBRyh(filename,location,whereput)
% 
% inputs are 
% filename: complete path and name of raw .rsk file
%  location: string describing where teh data were collected (used in name
%  creation)
% whereput is name of output directory where to save stuff
% netcdf_out = name of saved netcdf file including location
% 
% [optional]
% -has option to save to standard text file; commented out by default
% 
% [modifications]
% 20180530 
% -changed variable for name so that works with new rsk tools ( but this
% then failed to read teh temperature data... so use old tools instead)
% - commented out save to text file bc redundant
% - made it save nc file into processing directory
% 20211125
% -started using newer 3.5.3 versoin of rsk tools to enable reading of new
% format of files saved in ruskin.
% 20220921 yh
% added writing UTC offset to nc file
% added verbosity
%%
% % clear all
% % close all
% % clc;
% 
% % javaaddpath('/Users/Yasha/Java/netcdfAll-4.3.jar');
% % javaaddpath ('/Users/Yasha/Documents/MATLAB/m-files/mexcdf/snctools/classes');
% % addpath ('/Users/Yasha/Documents/MATLAB/m-files/mexcdf/mexnc');
% % addpath ( '/Users/Yasha/Documents/MATLAB/m-files/mexcdf/snctools');
% % addpath ( '/Users/Yasha/Documents/MATLAB/m-files/matlab_netCDF_OPeNDAP');
% % addpath ( '/Users/Yasha/Documents/MATLAB/m-files/m_map');
% % addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/COLORMAPS_BARS'));
% % addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/jlab'));
% % addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rotateXLabels'));
% % addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/subaxis'));
% % addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files'));
% % addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/subaxis'));
% % addpath ('/Users/Yasha/Documents/MATLAB/m-files/peakfinder:');
% 
% % addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rbr-rsktools'));
% %% settings
% % 
% % 
% % fname='/Volumes/Margs/RBRpressure_sensors/SN77825_Gracetown_Jan2016/Gracetown_boatramp_077825_20160424_2037.rsk'
% % fname='Gracetown_boatramp_077825_20160424_2037.rsk'
% % fname='/Volumes/Margs/RBRpressure_sensors/SN124028_Gracetown_20180513/124028_20180513_1011.rsk'
% % fname='/Users/Yasha/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/SN124029_Gracetown_offshore20180519_1928/SN124029_Gracetown_offshore20180519_1928.rsk'
% % loc='Gracetown_Offshore_v2';
% 
% 
% % fname='Gracetown_boatramp_077284_20160926_1002.rsk'
% % loc='Gracetown_boatramp';

%%
% organise filenames

ismac = ~isempty( regexp( whereput, '/', 'start' )); % figure out to use forwar dor back slashes
isPC = ~isempty( regexp( whereput, '\', 'start' ));
if isempty( regexp( whereput, '/', 'start' )) || isempty( regexp( whereput, '\', 'start' ));
    nopath=0;
end

if nopath>0
    filebase=filename(1:end-4);
elseif ismac>0
naminds=regexp( filename, '/', 'start' );
elseif isPC>0
naminds=regexp( filename, '\', 'start' );
end

fname=[filename];
loc=location;
filebase=[filename(naminds(length(naminds))+1:length(filename)-4) '_' loc];



%% read and plot 



tic

data1=RSKopen(fname)
% RSKplotthumbnail(data1) % plot subsample
RSKplotdownsample(data1)
disp('Reading complete rsk file...SLOW')
data1=RSKreaddata(data1); % read all data


if strcmpi('OFFSET_FROM_UTC',data1.parameterKeys(16).key)
    OFFSET_FROM_UTC=data1.parameterKeys(16).value
end

%subtract standard atm pressure or use data
% if exist('mslp_time','var')
%     corrected_with_mslp_obs='Y';
%     
%     it_mslp=data1.data.tstamp;
%     imslp=interp1(mslp_time,mslp,it_mslp);
%     
%     data1.data.values(:,1)=data1.data.values(:,1)-10.1325;
%     
% else
%     corrected_with_mslp_obs='N';
% data1.data.values(:,1)=data1.data.values(:,1)-10.1325;
data1.data.values(:,1)=data1.data.values(:,1);
% end


% 
% % plot subset of data
% figure;
% hold on
% plot(data1.thumbnailData.tstamp,data1.thumbnailData.values(:,1)-mean(data1.thumbnailData.values(:,1)),'r')
% % legend('huzza','ramp')
% datetick


% plot all data
figure;
hold on
plot(data1.data.tstamp,data1.data.values(:,1)-mean(data1.data.values(:,1)),'r') %if using new rsk tools
% plot(data1.data.tstamp,data1.datasets.values(:,1)-mean(data1.data.values(:,1)),'r') % old rsktools
datetick

%% subset to good data
% 

t1=data1.data.tstamp;
p1=data1.data.values(:,1);


if size(data1.data.values,2)==2

temp1=data1.data.values(:,2);

end



%% save text [optional]

% % % % % % ismac = ~isempty( regexp( whereput, '/', 'start' ))
% % % % % % isPC = ~isempty( regexp( whereput, '\', 'start' ))
% 
% if ismac>0
% fnc=[whereput '/' filebase  '.txt'];
% elseif isPC>0
% fnc=[whereput '\' filebase  '.txt'];  
% else
%     disp('No path given... saving text file to current directory')
%     fnc=[filebase  '.txt'];
% end
% % % % % % 
% % % % % % fn=[whereput '/' fname(1:end-4) '.txt'];
% % % % % % fn=[fname(1:end-4) '.txt'];
% % % % 
% % % % 
% % % % % fid = fopen(fn,'w');
% % % % % 
% % % % % if exist('temp1','var')
% % % % % disp('Temperature exists')
% % % % % fprintf(fid,[data1.instruments.model ', S/N ' num2str(data1.instruments.serialID)  ', Raw file: ' data1.deployments.name  ', Time in units of days since 2000-01-01, sea pressure (column 1)  in db is abs pressure -10.1325; temperature (column 2) in deg. C \n']);
% % % % % fprintf(fid, '%0.9f\t%0.5f\t%0.5f\n', [t1 p1 temp1]');
% % % % % fclose(fid)
% % % % % 
% % % % % else
% % % % %     
% % % % % fprintf(fid,[data1.instruments.model ', S/N ' num2str(data1.instruments.serialID)  ', Raw file: ' data1.deployments.name  ', Time in units of days since 2000-01-01, sea pressure (column 1)  in db is abs pressure -10.1325 \n']);
% % % % % fprintf(fid, '%0.9f\t%0.5f\n', [t1 p1 ]');
% % % % % fclose(fid)
% % % % %     
% % % % % end 
% % % % %     
% % % % % text_out=fn;
%% write to netcdf instead (need to fix to conform to normal style)
disp('Writing NETCDF file...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edit here for different files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fnc=[fname(1:end-4) '.nc'];
% ismac = ~isempty( regexp( whereput, '/', 'start' ))
% isPC = ~isempty( regexp( whereput, '\', 'start' ))

if ismac>0
fnc=[whereput '/' filebase  '.nc'];
elseif isPC>0
fnc=[whereput '\' filebase  '.nc'];  
else
    disp('No path given... saving nc file to current directory')
    fnc=[filebase  '.nc'];
end

p=p1;
% % t=t1; % this is matlab time
t=t1-datenum(2000,1,1); % this is days since 2000
if exist('temp1','var')
temp=temp1;
end
% t=t+datenum(2000,1,1); % this is back to serial time
% t=t-datenum(2000,1,1);% this is days since 2000
model=data1.instruments.model ;
serialID=num2str(data1.instruments.serialID) ;
Rawfilename=data1.deployments.name;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create file here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nccreate(fnc,'press','Dimensions',{'time',length(p)});
if exist('temp1','var')
nccreate(fnc,'temperature','Dimensions',{'time',length(p)});
end
ncwriteatt(fnc,'/','creation_time',datestr(now));
ncwriteatt(fnc,'/','Author','Yasha Hetzel');
ncwriteatt(fnc,'/','Location',loc);
ncwriteatt(fnc,'/','Data_period',[datestr(t1(1)) ' to ' datestr(t1(end))]);
if exist('OFFSET_FROM_UTC','var')
ncwriteatt(fnc,'/','Offset from UTC', OFFSET_FROM_UTC);
end
ncwriteatt(fnc,'/','Instrument_model',model);
ncwriteatt(fnc,'/','Serial_Number',serialID);
ncwriteatt(fnc,'/','Raw_filename',Rawfilename);

ncwrite(fnc,'press',p);
ncwriteatt(fnc,'press','units','Decibars');
ncwriteatt(fnc,'press','description','To get absolute pressure subtract 10.1325');
ncwriteatt(fnc,'press','description','Includes atm pressure. Must either subtract 10.1325 dB or use air pressure observations');

if exist('temp1','var')
ncwrite(fnc,'temperature',temp);
ncwriteatt(fnc,'temperature','units','Degrees C');
end

nccreate(fnc , 'time','dimensions',{'time',length(t)});
ncwriteatt(fnc,'time','units','Days since 2000-01-01T00:00:00');
if exist('OFFSET_FROM_UTC','var')
ncwriteatt(fnc,'time','Offset_from_UTC', OFFSET_FROM_UTC);
end
ncwrite(fnc,'time',t);

ncdisp(fnc);
toc

% save output filenames
netcdf_out=fnc;


disp(['NETCDF FILE SAVED HERE: ' netcdf_out])
end