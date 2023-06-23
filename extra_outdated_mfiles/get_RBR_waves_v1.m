
restoredefaultpath
addpath('/Users/Yasha/Documents/MATLAB/m-files')
addpath ( '/Users/Yasha/Documents/MATLAB/m-files/mexcdf/snctools');


% ncref='/Volumes/Margs/RBRpressure_sensors/SN77824_Gracetown_boatramp_July2016/SN77824_Gracetown_boatramp_July2016.nc'

ncref='/Volumes/Margs/RBRpressure_sensors/SN77824_Gracetown_boatramp_July2016/SN77825_Gracetown_Huzza_July2016.nc'

ts=nc_varget(ncref,'time')+datenum(2000,1,1);
p=nc_varget(ncref,'press');

clear Hs Ts Htime
numdays=floor(ts(end)-ts(1));
clear waves
co=0;
for hr=1:numdays*24
   co=co+1; 
   inds=[hr*7200-7199:hr*7200];
psub=p(inds);
tsub=ts(inds);
clear waves
waves=wavepar_yh(psub,2,1/30);

Hs(co)=waves.Hs;
Ts(co)=waves.Ts;
Htime(co)=tsub(3600);
end

figure; 
plot(Htime,Hs)
title('Hs')
datetick

figure; 
plot(Htime,Ts)
title('Ts')
datetick