
% yasha hetzel nov 2022
% read and calculate waves for chari from Cockburn Sound Stirling channel
% marker data from Fremantle Ports
% should run from within yasha-RBR_tools_v7 directory
% requirements: import_pressure_logfiles.m

% need to consider if




clear all;
clc;
close all;



restoredefaultpath % this just fixes if there are any conflicting functions


% % this should work if you run it from the yasha_RBR_tools directory
% addpath(genpath(pwd))
% addpath(genpath('rbr-rsktools'))

%%-------------------------------------------------------------------------%
%%                          Settings
%%-------------------------------------------------------------------------%

% if want to loop through files use this
fdir='/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/'
d=dir([fdir '*.log'])

for fi=1:length(d);
    filename=[d(fi).folder filesep d(fi).name]
%     
    
    
% filename='/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/2017-07-06T10-00.log';
% filename='/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/2018-08-06T15-00.log';
%%-------------------------------------------------------------------------%
%                       read data
%%-------------------------------------------------------------------------%
% [dat serialtime press] = import_pressure_logfiles(filename);

[dat ts p] = import_pressure_logfiles(filename);

p=p-10.1325; % remove mean atm pressure. It would be better if you had observations. I have them for this period but there are a few gaps so probably not best. If you are just calculating waves it doesnt matter.


%%-------------------------------------------------------------------------%
%            Wave calculation settings [for chari methods]
%             
%%-------------------------------------------------------------------------%
calculate_waves='chari_method' % [comment out if using ocn toolbox]
srate=2; % sample rate in hertz  [for chari method]
burst_duration= length(ts)/srate %2048 %3598 % 3600 %length(ts)/srate
% calculate_waves='ocn_toolbox'  [comment out if using chari method]

%%-------------------------------------------------------------------------%
%      Wave calculation settings [for oceanlyz toolbox]
%             default settings can be left 'as is'
%%-------------------------------------------------------------------------%
calculate_waves='ocn_toolbox'
InputType='waterlevel' %'pressure'              % Data input ['pressure'] or  'waterlevel'. If pressure, the signal attenuation with depth will be applied [recommended].
OutputType='wave' %'wave+waterlevel' %'wave+waterlevel'      % ['wave'], ['wave+waterlevel'];
AnalysisMethod='spectral'         % Wave calculation method. 'zerocross' or 'spectral'. use zerocross if you want Hs, but very similar to Hm0 and swell/sea are useful
fs=2;                             % instrument sample rate  in Hz
burst_duration= length(ts)/fs %3599 %3600 15*60 ;       % Interval of time window (in seconds) over which to calculate waves (must be <= number of samples in file)
site_depth=2;                      % approx depth of site (m) (seabed to surface)
heightfrombed=0 %site_depth-1 ;  % [Default=0]   % height of instrument above bed (m). Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'.If unknown assume = 0
dispout='no';                     % Display output of calculations at each burst ('yes or 'no')
Rho=1024;                         % Seawater density (Varies)
SeparateSeaSwell='yes'            % 'yes or 'no' (only works with spectral)
fmaxswell=1/8 %0.25;              % Maximum swell frequency (spectral only). 1/Period (ie. minimum period a swell can have). It should be between 0 and (fs/2). Only required if SeparateSeaSwell='yes' and AnalysisMethod='spectral'
fpminswell=1/25 %0.1;             % Minimum swell frequency (spectral only)
max_period= 25;                   % max period to use in wave calculations and filtering


%%-------------------------------------------------------------------------%
%                      Plotting settings
%%-------------------------------------------------------------------------%
plot_figs='N'  
save_matfig= 'N'
print_fig= 'N'
use_export_fig= 'N'
img_type='png'
fname   = 'test' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Use OceanLyz wave toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch calculate_waves
    case 'ocn_toolbox'

        n_bursts=floor(length(ts)./(burst_duration*fs));
        samples=burst_duration*fs;
        
        co=0;
        for i=1:n_bursts
            
            i1=co*samples+1;
            i2=i1+samples-1;
            burst_times(i)=mean(ts(i1:i2));
            co=co+1;
        end
        
        %add path for OCEANLYZ folder
        addpath(genpath('oceanlyz_2_0'));
        % cd('oceanlyz_2_0')
        
        %Create OCEANLYZ object
        clear ocn %Optional
        ocn=oceanlyz;
        
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
        ocn.heightfrombed=heightfrombed;   % %0.05;           %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
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
        
        ocn.wave
 
%  if you want to save in a loop must save each output here       
% %         wavedata(fi)=ocn.wave;
% % %         wavedata(fi).Eta=[];wavedata(fi).Burst_Data=[];
% %         wavetimes(fi)=burst_times;

end

%%-------------------------------------------------------------------------%
%% ------------------------------------------------------------------------%

calculate_waves='chari_method' % force to also do chari way
switch calculate_waves
    case 'chari_method'
%%-------------------------------------------------------------------------%        
   

% % % [wavedata hs1  hm1 S Pstat]=get_chari_wave_analysis(ts,p,srate);

[wavedata hs1  hm1 S Pn Pstat]=get_chari_wave_analysis(ts,p,srate,burst_duration);
Hs=std(p)*4; % Sig wave height

close
% display wave results
vars=  {'Tpeak' 'Tot_Var' 'Trend_Var' 'Inf_Var' 'Swell_Var' 'Wind_Var' 'Noise_Var' ...
	    'Trend_Per' 'Inf_Per' 'Swell_Per' 'Wind_Per' 'Noise_Per' ...
	   'Hrms' 'Trend_ht' 'Inf_ht' 'Swell_ht' 'Wind_ht' 'Noise_ht'};
   disp('-------------------------------------')
   disp(d(fi).name)
   disp(['Time: ' datestr(mean(ts))])
   disp('-------------------------------------')
for k=1:length(Pstat); disp ([num2str(k) '. ' vars{k} ' ' num2str(Pstat(k))]);end
disp ([num2str(k+1) '. Hs ' num2str(Hs)])

end
      
%% ------------------------------------------------------------------------%
% ----------------------------Plots-------------------------------------- %
%%-------------------------------------------------------------------------%

%         
% switch plot_figs
%     case 'Y'
%         % make sure time is correct length
%         burst_times=burst_times(1:size(ocn.wave.Burst_Data,1));
%         
%         xlims=[ts(1) ts(end)]; % all
%         xlims=[begin finish];% good
%         xt=[xlims(1):xlims(2)];
%         yt=[0:2:30];
%         %
%         % sp_corrected1=ocn.wave;
%         % sep_swell= ocn.SeparateSeaSwell;
%         switch AnalysisMethod
%             case 'spectral'
%                 
%                 switch SeparateSeaSwell
%                     case 'yes'
%                         %         figure;
%                         %         hold on;
%                         %         plot(ocn.wave.Hm0swell);
%                         %         plot(ocn.wave.Hm0sea);
%                         %         plot(ocn.wave.Hm0,'k');
%                         %         legend('swell','sea','Hm0');
%                         %
%                         %         figure;
%                         %         hold on;
%                         %
%                         %         h1=plot(ocn.wave.Tp);
%                         %         h2=plot(ocn.wave.Tpswell);
%                         %         h3=plot(ocn.wave.Tpsea);
%                         %         legend('Tp','swell','sea');
%                         %         ylabel('Period (s)')
%                         %
%                         %         h1.Color=[0 0 0];
%                         %         h2.Color=[ 0    0.4470    0.7410];
%                         %         h3.Color=[0.8500    0.3250    0.0980];
%                         
%                         f=figure;
%                         f.Position=[66 207 1600 450];
%                         f.Color='w';
%                         % set(gcf,'Position',[200 100 700 300])
%                         % set(gcf,'PaperPositionMode','auto')
%                         plot(burst_times,ocn.wave.Hm0swell)
%                         hold on
%                         plot(burst_times,ocn.wave.Hm0sea)
%                         plot(burst_times,ocn.wave.Hm0,'k')
%                         legend('Hm0swell','Hm0sea','Hm0')
%                         ylabel('Wave height (m)')
%                         set(gca,'xlim',xlims,'xtick',xt)
%                         set(gca,'ylim',[0 ceil(max(ocn.wave.Hm0))])
%                         datetick('x','mm/dd','keeplimits','keepticks')
%                         rotateXLabels( gca(), 90 )
%                         xlabel(datestr(mean(burst_times),'YYYY'))
%                         title([instrument ' ' location ],'interpreter','none')
%                         fname=[whereput '/' instrument '_' location '_WAVE_HEIGHT_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
%                         % print('-dpng','-r200', fname);
%                         % print('-dpdf', '-painters', fname);
%                         %         export_fig(fname,'-pdf','-png','-transparent')
%                         savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
%                         
%                         f=figure;
%                         f.Position=[200 100 1600 450];
%                         f.Color='w';
%                         h1=plot(burst_times,ocn.wave.Tp)
%                         hold on
%                         h2=plot(burst_times,ocn.wave.Tpswell)
%                         h3=plot(burst_times,ocn.wave.Tpsea)
%                         set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
%                         legend('Tp','swell','sea');
%                         ylabel('Period (s)')
%                         title([instrument ' ' location ],'interpreter','none')
%                         datetick('x','mm/dd','keeplimits','keepticks')
%                         %         xlabel(datestr(mean(burst_times),'YYYY'))
%                         if exist('Offset_from_UTC','var')
%                             xlabel([datestr(mean(burst_times),'YYYY') ' (UTC+' Offset_from_UTC ')'])
%                         else
%                             xlabel(datestr(mean(burst_times),'YYYY'))
%                         end
%                         rotateXLabels( gca(), 90 )
%                         h1.Color=[0 0 0];
%                         h2.Color=[ 0    0.4470    0.7410];
%                         h3.Color=[0.8500    0.3250    0.0980];
%                         
%                         fname=[whereput '/' instrument '_' location '_WAVE_PERIOD_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
%                         % print('-dpng','-r200', fname);
%                         % print('-dpdf', '-painters', fname);
%                         %         export_fig(fname,'-pdf','-png','-transparent')
%                         savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
%                         
%                         
%                     case 'no'
%                         
%                         %                     figure;
%                         %                     hold on;
%                         %                     plot(ocn.wave.Hm0);
%                         %                     legend('Hm0');
%                         %
%                         %                     figure;
%                         %                     hold on;
%                         %                     plot(ocn.wave.Tp);
%                         %                     legend('Tp');
%                         %                     ylabel('Period (s)')
%                         
%                         f=figure;
%                         f.Position=[66 207 1600 450];
%                         f.Color='w';
%                         plot(burst_times,ocn.wave.Hm0,'k')
%                         legend('Hm0')
%                         ylabel('Wave height (m)')
%                         set(gca,'xlim',xlims,'xtick',xt)
%                         set(gca,'ylim',[0 ceil(max(ocn.wave.Hm0))])
%                         datetick('x','mm/dd','keeplimits','keepticks')
%                         rotateXLabels( gca(), 90 )
%                         xlabel(datestr(mean(burst_times),'YYYY'))
%                         title([instrument ' ' location ],'interpreter','none')
%                         fname=[whereput '/' instrument '_' location '_WAVE_HEIGHT_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
%                         % print('-dpng','-r200', fname);
%                         % print('-dpdf', '-painters', fname);
%                         %         export_fig(fname,'-pdf','-png','-transparent')
%                         savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
%                         
%                         f=figure;
%                         f.Position=[200 100 1600 450];
%                         f.Color='w';
%                         h1=plot(burst_times,ocn.wave.Tp)
%                         set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
%                         legend('Tp','swell','sea');
%                         ylabel('Period (s)')
%                         title([instrument ' ' location ],'interpreter','none')
%                         datetick('x','mm/dd','keeplimits','keepticks')
%                         %         xlabel(datestr(mean(burst_times),'YYYY'))
%                         if exist('Offset_from_UTC','var')
%                             xlabel([datestr(mean(burst_times),'YYYY') ' (UTC+' Offset_from_UTC ')'])
%                         else
%                             xlabel(datestr(mean(burst_times),'YYYY'))
%                         end
%                         rotateXLabels( gca(), 90 )
%                         h1.Color=[0 0 0];
%                         %                 h2.Color=[ 0    0.4470    0.7410];
%                         %                 h3.Color=[0.8500    0.3250    0.0980];
%                         
%                         fname=[whereput '/' instrument '_' location '_WAVE_PERIOD_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
%                         % print('-dpng','-r200', fname);
%                         % print('-dpdf', '-painters', fname);
%                         %         export_fig(fname,'-pdf','-png','-transparent')
%                         savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
%                         
%                         
%                 end
%                 
%                 
%             case 'zerocross'
%                 %              Hs: [100×1 double]
%                 %              Hz: [100×1 double]
%                 %              Tz: [100×1 double]
%                 %              Ts: [100×1 double]
%                 % Hmax
%                 
%                 f=figure;
%                 f.Position=[66 207 1600 450];
%                 f.Color='w';
%                 % set(gcf,'Position',[200 100 700 300])
%                 % set(gcf,'PaperPositionMode','auto')
%                 plot(burst_times,ocn.wave.Hs)
%                 hold on
%                 plot(burst_times,ocn.wave.Hz)
%                 plot(burst_times,ocn.wave.Hmax)
%                 legend('Hs','Hz','Hmax')
%                 ylabel('Wave height (m)')
%                 set(gca,'xlim',xlims,'xtick',xt)
%                 set(gca,'ylim',[0 ceil(max(ocn.wave.Hmax))])
%                 datetick('x','mm/dd','keeplimits','keepticks')
%                 rotateXLabels( gca(), 90 )
%                 xlabel(datestr(mean(burst_times),'YYYY'))
%                 title([instrument ' ' location ],'interpreter','none')
%                 fname=[whereput '/' instrument '_' location '_WAVE_HEIGHT_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
%                 % print('-dpng','-r200', fname);
%                 % print('-dpdf', '-painters', fname);
%                 %         export_fig(fname,'-pdf','-png','-transparent')
%                 savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
%                 
%                 f=figure;
%                 f.Position=[200 500 1600 450];
%                 f.Color='w';
%                 plot(burst_times,ocn.wave.Ts)
%                 hold on
%                 plot(burst_times,ocn.wave.Tz)
%                 set(gca,'xlim',xlims,'xtick',xt,'ytick',yt)
%                 legend('Ts','Tz')
%                 ylabel('Period (s)')
%                 title([instrument ' ' location ],'interpreter','none')
%                 datetick('x','mm/dd','keeplimits','keepticks')
%                 %         xlabel(datestr(mean(burst_times),'YYYY'))
%                 if exist('Offset_from_UTC','var')
%                     xlabel([datestr(mean(burst_times),'YYYY') ' (UTC+' Offset_from_UTC ')'])
%                 else
%                     xlabel(datestr(mean(burst_times),'YYYY'))
%                 end
%                 rotateXLabels( gca(), 90 )
%                 fname=[whereput '/' instrument '_' location '_WAVE_PERIOD_' datestr(xlims(1),'yyyymmdd') '-' datestr(xlims(2),'yyyymmdd')];
%                 % print('-dpng','-r200', fname);
%                 % print('-dpdf', '-painters', fname);
%                 %         export_fig(fname,'-pdf','-png','-transparent')
%                 savefig2_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
%                 
%         end
%         
% end     
%         


end % loop through files
