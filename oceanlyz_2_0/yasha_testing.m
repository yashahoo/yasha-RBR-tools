            

% yasha  testing oceanlyz

%% working loop no pressure attentuation correction ( must input corrected pressure or  water level data)

ncref='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/208787_20220911_1406_BusseltonJetty341_end_processed/208787_20220911_1406_BusseltonJetty341_end.nc'
press=ncread(ncref,'press');
time=ncread(ncref,'time')+datenum(2000,1,1);


burst_duration=3600; % secs
fs=2; % Hz

n_bursts=floor(length(time)./(burst_duration*fs));



co=0
for i=1:n_bursts
i1=co*burst_duration+1;
i2=i1+burst_duration*fs-1;
    tmp=press(i1:i2);
[Hs(i),Hz(i),Tz(i),Ts(i),H,T]=WaveZerocrossingFun(tmp,fs,burst_duration,'off');
Hmax(i)=max(H);
wTime(i)=mean(time(i1:i2));
co=co+1

end



figure; 
plot(wTime,Hs);
hold on
plot(wTime,Hz);
plot(wTime,Hmax);
datetick


%% RUN  AS FUNCTION
ncref='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/208787_20220911_1406_BusseltonJetty341_end_processed/208787_20220911_1406_BusseltonJetty341_end.nc'
press=ncread(ncref,'press');
time=ncread(ncref,'time')+datenum(2000,1,1);
addpath(genpath('oceanlyz_2_0'))
[Hs,Hz,Tz,Ts,Hmax,Tmax,wTime] = get_all_WaveZerocrossingFun_fn(time(1:3600*2*100),press(1:3600*2*100),2,3600);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try using whole package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run with zero cross (waterlevel)
    
    
% ncref='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/208787_20220911_1406_BusseltonJetty341_end_processed/208787_20220911_1406_BusseltonJetty341_end.nc'
ncref='/Volumes/Margs_Clone/RBRpressure_sensors/SN124028_Gracetown_offshore_20180819_2028/124028_20180819_2028_processed/124028_20180819_2028.nc'
press=ncread(ncref,'press');
time=ncread(ncref,'time');%+datenum(2000,1,1);

burst_duration=3600;
fs=2;
n_bursts=floor(length(time)./(burst_duration*fs));

% n_bursts=100;
    
    
        %Assume OCEANLYZ files are in 'C:\oceanlyz_matlab' folder
    cd('/Users/00068592/Documents/MATLAB/m-files/yashatools/yasha_RBR_tools_v5/oceanlyz_2_0')

    %Create OCEANLYZ object
    %clear ocn %Optional
    ocn=oceanlyz;
    
    %Read data
    %Assume data file is named 'waterpressure_5burst.csv' and is stored in 'C:\oceanlyz_matlab\Sample_Data'
    current_folder=pwd;                  %Current (OCEANLYZ) path
    cd('Sample_Data/') %Change current path to Sample_Data folder
    water_pressure=press;%+10.13;%*10000; %Load data
    cd(current_folder)                   %Change current path to OCEANLYZ folder
%     water_pressure=water_pressure(1:n_bursts*3600*fs);
    
    %Input parameters
    ocn.data= water_pressure;
    ocn.InputType='waterlevel';
    ocn.OutputType='wave' %'wave+waterlevel';
    ocn.AnalysisMethod='zerocross';
    ocn.n_burst=n_bursts;
    ocn.burst_duration=3600;
    ocn.fs=2;
    ocn.fmin=0.05;                    %Only required if ocn.AnalysisMethod='spectral'
    ocn.fmax=ocn.fs/2;                %Only required if ocn.AnalysisMethod='spectral'
    ocn.fmaxpcorrCalcMethod='auto';   %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.Kpafterfmaxpcorr='constant';  %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.fminpcorr=0.15;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.fmaxpcorr=0.55;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.heightfrombed=0 %0.05;           %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.dispout='no';
    ocn.Rho=1024;                     %Seawater density (Varies)

    %Run OCEANLYZ
    ocn.runoceanlyz()

    %Plot peak wave period (Tp)
    figure;
    plot(ocn.wave.Hs)
    
zc_no_correction=ocn.wave;
%% run with zero cross (pressure and correction)
    
    
% ncref='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/208787_20220911_1406_BusseltonJetty341_end_processed/208787_20220911_1406_BusseltonJetty341_end.nc'
ncref='/Volumes/Margs_Clone/RBRpressure_sensors/SN124028_Gracetown_offshore_20180819_2028/124028_20180819_2028_processed/124028_20180819_2028.nc'
press=ncread(ncref,'press');
time=ncread(ncref,'time')+datenum(2000,1,1);

burst_duration=3600;
fs=2;
n_bursts=floor(length(time)./(burst_duration*fs));
n_bursts=100;
% n_bursts=50;
    
    
        %Assume OCEANLYZ files are in 'C:\oceanlyz_matlab' folder
    cd('/Users/00068592/Documents/MATLAB/m-files/yashatools/yasha_RBR_tools_v5/oceanlyz_2_0')

    %Create OCEANLYZ object
    %clear ocn %Optional
    ocn=oceanlyz;
    
    %Read data
    %Assume data file is named 'waterpressure_5burst.csv' and is stored in 'C:\oceanlyz_matlab\Sample_Data'
    current_folder=pwd;                  %Current (OCEANLYZ) path
    water_pressure=10000*(press);%+10.13);%*10000; %Load data convert to n/m2 (+3010 to test for depth)
    cd(current_folder)                   %Change current path to OCEANLYZ folder
    water_pressure=water_pressure(1:n_bursts*3600*fs); % not needed
    
    %Input parameters
    ocn.data= water_pressure;
    ocn.InputType='pressure' %'waterlevel';
    ocn.OutputType='wave+waterlevel';
    ocn.AnalysisMethod='zerocross';
    ocn.n_burst=n_bursts;
    ocn.burst_duration=3600;
    ocn.fs=fs;
%     ocn.fmin=0.05;                    %Only required if ocn.AnalysisMethod='spectral'
%     ocn.fmax=ocn.fs/2;                %Only required if ocn.AnalysisMethod='spectral'
%     ocn.fmaxpcorrCalcMethod='auto';   %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
%     ocn.Kpafterfmaxpcorr='constant';  %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
%     ocn.fminpcorr=0.15;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
%     ocn.fmaxpcorr=0.55;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.heightfrombed=0 %6.5%0.05;           %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.dispout='no';
    ocn.Rho=1024;                     %Seawater density (Varies)

    %Run OCEANLYZ
    ocn.runoceanlyz()

    %Plot peak wave period (Tp)
    figure;
    plot(ocn.wave.Hs)

    
     
zp_corrected3=ocn;
%% spectral approach
% ncref='/Volumes/Margs_Clone/RBRpressure_sensors/Busselton_Jetty_202208/208787_20220911_1406_BusseltonJetty341_end_processed/208787_20220911_1406_BusseltonJetty341_end.nc'
ncref='/Volumes/Margs_Clone/RBRpressure_sensors/SN124028_Gracetown_offshore_20180819_2028/124028_20180819_2028_processed/124028_20180819_2028.nc'

press=ncread(ncref,'press');
time=ncread(ncref,'time');%+datenum(2000,1,1);

burst_duration=3600;
fs=2;
n_bursts=floor(length(time)./(burst_duration*fs));
btime=[time(1):burst_duration/3600/24:time(n_bursts*burst_duration)];
% n_bursts=100 %

%Change current working directory to OCEANLYZ folder
%Assume OCEANLYZ files are in 'C:\oceanlyz_matlab' folder
cd('/Users/00068592/Documents/MATLAB/m-files/yashatools/yasha_RBR_tools_v5/oceanlyz_2_0')

%Create OCEANLYZ object
%clear ocn %Optional
ocn=oceanlyz;

%Read data
current_folder=pwd;                  %Current (OCEANLYZ) path
water_pressure=press*10000; %Load data
cd(current_folder)                   %Change current path to OCEANLYZ folder
water_pressure=water_pressure(1:n_bursts*3600*fs);




%Input parameters
ocn.data= water_pressure;
ocn.InputType='pressure' %'waterlevel';
ocn.OutputType='wave' %'wave+waterlevel';
ocn.AnalysisMethod='spectral';
ocn.n_burst=n_bursts;
ocn.burst_duration=3600;
ocn.fs=2;
%     ocn.fmin=1/20
%     ocn.fmax=ocn.fs/2;  %0.1250
ocn.fmin=0.05;                    %Only required if ocn.AnalysisMethod='spectral'
ocn.fmax=ocn.fs/2;                %Only required if ocn.AnalysisMethod='spectral'
ocn.fmaxpcorrCalcMethod='auto';   %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
ocn.Kpafterfmaxpcorr='constant';  %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
ocn.fminpcorr=0.15;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
ocn.fmaxpcorr=0.55;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
ocn.heightfrombed=0 %0.05;           %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
ocn.dispout='no';
ocn.Rho=1024;                     %Seawater density (Varies)

ocn.SeparateSeaSwell='yes'
%maximum swell frequency
ocn.fmaxswell=1/8 %0.25;
%                                 Maximum frequency that swell can have (It is about 0.2 in Gulf of Mexico) in (Hz)
%                                     It should be between 0 and (fs/2)
%                                     Only required if SeparateSeaSwell='yes' and AnalysisMethod='spectral'

%Minimum swell frequency
ocn.fpminswell=1/25 %0.1;
%                                 Minimum frequency that swell can have (it is used for Tpswell calculation) in (Hz)
%                                     It should be between 0 and (fs/2)
%                                     Only required if SeparateSeaSwell='yes' and AnalysisMethod='spectral'


%Run OCEANLYZ
ocn.runoceanlyz()

%   %Plots
sp_corrected1=ocn.wave;
sep_swell= ocn.SeparateSeaSwell;

switch sep_swell
    case 'yes'
        figure;
        hold on;
        plot(ocn.wave.Hm0swell);
        plot(ocn.wave.Hm0sea);
        plot(ocn.wave.Hm0,'k');
        legend('swell','sea','Hm0');
        
        figure;
        hold on;
        
        h1=plot(ocn.wave.Tp);
        h2=plot(ocn.wave.Tpswell);
        h3=plot(ocn.wave.Tpsea);
        legend('Tp','swell','sea');
        ylabel('Period (s)')
        
        h1.Color=[0 0 0];
        h2.Color=[ 0    0.4470    0.7410];
        h3.Color=[0.8500    0.3250    0.0980];
        
    case 'no'
        if strcmpi(ocn.AnalysisMethod, 'spectral')
            figure;
            hold on;
            plot(ocn.wave.Hm0);
            legend('Hm0');
            
            figure;
            hold on;
            plot(ocn.wave.Tp);
            legend('Tp');
            ylabel('Period (s)')
            
        else
            figure;
            hold on;
            plot(ocn.wave.Hs);
            legend('Hs');
            
            figure;
            hold on;
            plot(ocn.wave.Ts);
            legend('Ts');
            ylabel('Period (s)')
        end
end


   