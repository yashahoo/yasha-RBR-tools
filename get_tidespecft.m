function [logf logj timef]=get_tidespecft(time,data,numpoints,overlap,samplehertz);
% 
%Program to produce time-frequency plots from sea level data
% [logf logj]=get_tidespecft(time,data,numpoints,overlap,samplehertz);
% edited yh to make function
% close all
%clear

% load jurienSj4.mat   % load data file
% clear Sxx
% load R124024_KB_dp.mat
% get_tidespecft(ts,p,4096,1024,2);
% time=ts;
% data=p;
% numpoints=4096
% overlap=1024
% samplehertz=2


% dataq=data1(:,9);                            % assign pressure to variable data
dataq=data;
time1=time;

isectf=floor((length(dataq))/(overlap))-3;

% isectf=1000%
% isectf=12090;                                % number of sections to analyse
% %isectf=5644     
JJc=ones(59,isectf);

samplehertz=2.0; 
sampleperiod=1/samplehertz;

% numpoints=4096;
spoint=1;
% while epoint<length(dataq)
for i=1:isectf
    
    
epoint = spoint + numpoints - 1;

if epoint>length(dataq)
    disp(['Exiting loop at index: ' num2str(epoint)]);
    break
end

data1=dataq(spoint:epoint,1);
timef(i)=mean(time1(spoint:epoint,1));  % time is in matlab format
disp(['Analyzing ', num2str(numpoints), ' points from ', num2str(spoint), ' to ', num2str(epoint)])

% to produce spectral densities and frequencies

       out=oppsd2(sampleperiod,data1,1,1,64,1.095,1,1);
JJc(1:59,i)=out(1:59,2);
spoint = spoint+overlap;                       % create over lap
end
% end


f=out(:,1);
logf=log(f);
logj=log(JJc);

close
%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % To format the plot - Sxx vs f (log-log)
% clf reset
% %loglog(f,JJc)
% logf=log(f);
% logj=log(JJc);
% %xlim([735600 735965]) 
% %xlim([1 3960]) 
% pcolor(timef,logf,logj)
% shading('interp')
% datetick('x',19)
% hold on
% 
% %axis([735600 735965  -12 -6])
% % to add a vertical line at a period of 24 hours (f24)
% ymiddle=0.20;
% xlim=get(gca,'xlim');
% 
% f24=log(1/(12*60*60));
% plot(xlim,[f24 f24],'k-')
% text(timef(1000),f24+0.25,' 12h','fontsize',[12])
% 
% f24=log(1/(24*60*60));
% plot(xlim,[f24 f24],'k-')
% text(timef(1000),f24+0.25,' 24h','fontsize',[12])
% 
% f24=log(1/(6*60*60));
% %plot(xlim,[f24 f24],'k-')
% text(timef(1000),f24+0.25,'6hr','fontsize',[12])
% 
% f24=log(1/4.5);
% plot(xlim,[f24 f24],'k-')
% text(timef(1000),f24+0.25,' 4.5s','fontsize',[12])
% 
% f24=log(1/8);
% plot(xlim,[f24 f24],'k-')
% text(timef(1000),f24+0.25,' 8s','fontsize',[12])
% 
% f24=log(1/15);
% plot(xlim,[f24 f24],'k-')
% text(timef(1000),f24+0.25,' 15s','fontsize',[12])
% 
% f24=log(1/30);
% plot(xlim,[f24 f24],'k-')
% text(timef(1000),f24+0.25,' 30s','fontsize',[12])
% 
% % f24=log(1/130);
% % plot(xlim,[f24 f24],'k-')
% % text(timef(1000),f24+0.25,'130s','fontsize',[12])
% 
% % f24=log(1/200);
% % plot(xlim,[f24 f24],'k-')
% % text(timef(1000),f24+0.25,'200s','fontsize',[12])
% 
% 
% f24=log(1/300);
% plot(xlim,[f24 f24],'k-')
% text(timef(1000),f24+0.25,' 300s','fontsize',[12])
% xlabel('Time')
% ylabel('Frequency(Hz)')

