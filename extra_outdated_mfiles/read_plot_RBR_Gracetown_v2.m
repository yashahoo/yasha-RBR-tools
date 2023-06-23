% read rbr data with ruskin toolkit

%%
clear all
close all
clc;
javaaddpath('/Users/Yasha/Java/netcdfAll-4.3.jar');
javaaddpath ('/Users/Yasha/Documents/MATLAB/m-files/mexcdf/snctools/classes');
addpath ('/Users/Yasha/Documents/MATLAB/m-files/mexcdf/mexnc');
addpath ( '/Users/Yasha/Documents/MATLAB/m-files/mexcdf/snctools');
addpath ( '/Users/Yasha/Documents/MATLAB/m-files/matlab_netCDF_OPeNDAP');
addpath ( '/Users/Yasha/Documents/MATLAB/m-files/m_map');
addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/COLORMAPS_BARS'));
addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/jlab'));
addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rotateXLabels'));
addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/subaxis'));
addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files'));
addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/subaxis'));
addpath ('/Users/Yasha/Documents/MATLAB/m-files/peakfinder:');

addpath(genpath('/Users/Yasha/Documents/MATLAB/m-files/rbr-rsktools'));
%% settings
% 
% 
% fname='/Volumes/Margs/RBRpressure_sensors/SN77825_Gracetown_Jan2016/Gracetown_boatramp_077825_20160424_2037.rsk'
% fname='Gracetown_boatramp_077825_20160424_2037.rsk'
% fname='/Volumes/Margs/RBRpressure_sensors/SN124028_Gracetown_20180513/124028_20180513_1011.rsk'
fname='/Users/Yasha/Documents/RESEARCH/DATA/MEASURED/RBRpressure_sensors/SN124029_Gracetown_offshore20180519_1928/SN124029_Gracetown_offshore20180519_1928.rsk'
loc='Gracetown_Offshore_v2';


% fname='Gracetown_boatramp_077284_20160926_1002.rsk'
% loc='Gracetown_boatramp';

%% read and plot 


tic

data1=RSKopen(fname)
RSKplotthumbnail(data1) % plot subsample
data1=RSKreaddata(data1); % read all data
%subtract standard atm pressure
data1.data.values(:,1)=data1.data.values(:,1)-10.1325;



% plot subset of data
figure;
hold on
plot(data1.thumbnailData.tstamp,data1.thumbnailData.values(:,1)-mean(data1.thumbnailData.values(:,1)),'r')
% legend('huzza','ramp')
datetick


% plot all data
figure;
hold on
plot(data1.data.tstamp,data1.data.values(:,1)-mean(data1.data.values(:,1)),'r')
datetick

%% subset to good data
% 
t1=data1.data.tstamp;
p1=data1.data.values(:,1);


if size(data1.data.values,2)==2

temp1=data1.data.values(:,2);

end



%% save text
fn=[fname(1:end-4) '.txt'];
fid = fopen(fn,'w');

if exist('temp1','var')
disp('Temperature exists')
fprintf(fid,[data1.instruments.model ', S/N ' num2str(data1.instruments.serialID)  ', Raw file: ' data1.datasets.name  ', Time in units of days since 2000-01-01, sea pressure (column 1)  in db is abs pressure -10.1325; temperature (column 2) in deg. C \n']);
fprintf(fid, '%0.9f\t%0.5f\t%0.5f\n', [t1 p1 temp1]');
fclose(fid)

else
    
fprintf(fid,[data1.instruments.model ', S/N ' num2str(data1.instruments.serialID)  ', Raw file: ' data1.datasets.name  ', Time in units of days since 2000-01-01, sea pressure (column 1)  in db is abs pressure -10.1325 \n']);
fprintf(fid, '%0.9f\t%0.5f\n', [t1 p1 ]');
fclose(fid)
    
end 
    

%% write to netcdf instead (need to fix to conform to normal style)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edit here for different files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnc=[fname(1:end-4) '.nc'];
p=p1;
t=t1; % this is days since 2000
if exist('temp1','var')
temp=temp1;
end
% t=t+datenum(2000,1,1); % this is back to serial time
% t=t-datenum(2000,1,1);% this is days since 2000
model=data1.instruments.model ;
serialID=num2str(data1.instruments.serialID) ;
Rawfilename=data1.datasets.name;


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
ncwriteatt(fnc,'/','Data_period',[datestr(t(1)) ' to ' datestr(t(end))]);
ncwriteatt(fnc,'/','Instrument_model',model);
ncwriteatt(fnc,'/','Serial_Number',serialID);
ncwriteatt(fnc,'/','Raw_filename',Rawfilename);

ncwrite(fnc,'press',p);
ncwriteatt(fnc,'press','units','Decibars');
ncwriteatt(fnc,'press','description','to get absolute pressure subtract 10.1325');

if exist('temp1','var')
ncwrite(fnc,'temperature',temp);
ncwriteatt(fnc,'temperature','units','Degrees C');
end

nccreate(fnc , 'time','dimensions',{'time',length(t)});
ncwriteatt(fnc,'time','units','Days since 2000-01-01T00:00:00Z');
ncwrite(fnc,'time',t);

ncdisp(fnc);
toc