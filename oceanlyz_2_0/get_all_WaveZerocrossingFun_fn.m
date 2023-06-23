function [Hs,Hz,Tz,Ts,Hmax,Tmax,wTime] = get_all_WaveZerocrossingFun_fn(time,press,fs,burst_duration);

% [Hs,Hz,Tz,Ts,Hmax,wTime] = get_all_WaveZerocrossingFun_fn(time,press,fs,burst_duration);
% press = continuous timeseries of pressure or sealevel data
% fs = frequency in Hz
% burst_duration = processing interval for wave calculations in seconds (eg. 1800 or 3600)
% yasha hetzel 20221014

% burst_duration=3600; % secs
% fs=2; % Hz

n_bursts=floor(length(time)./(burst_duration*fs));



co=0; Hs=[];Hz=[];Tz=[];Ts=[];Hmax=[];wTime=[];Tmax=[];
for i=1:n_bursts
    i1=co*(burst_duration*fs)+1;
    i2=i1+(burst_duration*fs)-1;
    tmp=press(i1:i2);
    [Hs(i),Hz(i),Tz(i),Ts(i),H,T]=WaveZerocrossingFun(tmp,fs,burst_duration,'off');
    Hmax(i)=max(H);
    Tmax(i)=max(T);
    wTime(i)=mean(time(i1:i2));
    co=co+1
    
end



f=figure;
f.Position=[400 300 1000 700];
subplot(2,1,1)
plot(wTime,Hs);
hold on
plot(wTime,Hz);
plot(wTime,Hmax);
legend('Hs','Hz','Hmax')
datetick
set(gca,'fontsize',14)
ylabel('Wave height (m)')

subplot(2,1,2)
plot(wTime,Ts);
hold on
plot(wTime,Tz);
% plot(wTime,Tmax);
legend('Ts','Tz','Tmax')
datetick
set(gca,'fontsize',14)
ylabel('Period (s)')







end