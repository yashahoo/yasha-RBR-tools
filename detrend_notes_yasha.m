
dat=input(1:100000);
tt=ts(1:100000);

p = polyfit(1:length(dat),dat,1);
f = polyval(p,1:length(dat)); 

dt=detrend(dat);

% figure; hold on;
% plot(1:length(dat),dat)
% plot(tt,f)
% plot(tt,dat-f)
% plot(tt(1:50000),dt(1:50000))
% plot(tt(1:50000),dt(1:50000)+f(1:50000))
% legend('raw','trend','detrended with polyval','detrend function','fixed detrended')

figure; hold on;
plot(dat)
plot(f)
plot(dat-f)
plot(dt(1:50000))
plot(dt(1:50000)+f(1:50000))
legend('raw','trend','detrended with polyval','detrend function','fixed detrended')


%%

% yh edits to PcorZerocrossingFun_yh.m

% [input1, Tr]=detrend(input,'linear'); % need newer matlab
pf = polyfit(1:length(input),input,1); % yh
Tr = polyval(pf,1:length(input)); %yh
if isrow(Tr) & iscolumn(input)
    Tr=Tr';
elseif iscolumn(Tr) & isrow(input)
    Tr=Tr';
end

input1 = input - Tr; %yh


     
%  [Eta Tr]=PcorZerocrossingFun_yh(input,fs,duration,h,heightfrombed,dispout);
 figure; hold on;
plot(input)
plot(Tr)
plot(Eta)
plot(Eta(1:length(Eta)/2)+Tr(1:length(Eta)/2))
legend('raw','trend','detrended with polyval','fixed detrended')       

        


