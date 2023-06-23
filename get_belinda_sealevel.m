% yh 2023-03-17
% get belinda Fremantle sea level means

file='/Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/Belinda/penguin_paper/sealevel/FSL_1984-2019.xlsx'

B = xlsread(file, 'FSL_1984-2019');

time=datenum(B(:,1), B(:,2), B(:,3));
sl=B(:,4);
figure; plot(time,sl)
datetick

ind=find(time>datenum(1985,12,31) & time<datenum(2020,1,1));

t=time(ind);
z=sl(ind);
figure; plot(t,z)
datetick


mean(z,'omitnan')



zd=detrend(z,'omitnan');
figure; plot(t,zd)
datetick



idx = isnan(z);
len=1:length(z);
p = polyfit(len(~idx),z(~idx),1);
f = polyval(p,1:length(p)); 


figure; hold on;
plot(z,'k')
plot(f,'b')
plot(z-f,'r')
plot(zd,'g')
plot(zd+f,'c')
legend('raw','trend','detrended with polyval','detrend function','fixed detrended');


%%

[yy mm dd]=datevec(t);

mon=unique(mm);
yrs=unique(yy);
yco=0; mco=0;zm=[];zdm=[];zy=[];zdy=[];
for yr=yrs';
    yco=yco+1;
    disp(num2str(yr))
    yi=find(yy==yr);
    
    yt(yco)=datenum(yr,7,1);
    zy(yco)=mean(z(yi),'omitnan');
    zdy(yco)=mean(zd(yi),'omitnan');
    fy(yco)=mean(f(yi),'omitnan');
    
    for mo=mon';
        mco=mco+1;
        mi=find(yy==yr & mm==mo);
        zm(mco)=mean(z(mi),'omitnan');
        zdm(mco)=mean(zd(mi),'omitnan');
        mt(mco)=datenum(yr,mo,15);
        fm(mco)=mean(f(mi),'omitnan');
        
    end
end


% 3 month moving average
z3m  = movmean(zm,3,'omitnan');
zd3m = movmean(zdm,3,'omitnan');
%% plot it


 f=figure; 
 f.Position=[200 200 1600,900];
 f.Color='w';
 hold on; 

 plot(yt,zy);
 plot(mt,z3m);
 plot(yt,fy,'r');
 
 plot(yt,zdy);
 plot(mt,zd3m);
 legend('annual mean','3mo moving average','trend','detrended annual','detrended 3mo moving average');
  ylabel('Fremantle Sealevel (mm)')
  set(gca,'xlim',[yt(1) mt(end)]);
  datetick('x','yyyy','keeplimits')
  set(gca,'fontsize',14)
  grid on
  box on





 f=figure; 
 f.Position=[200 200 1600,900];
 f.Color='w';
 hold on; 
plot(t,z,'color',[.8 .8 .8]);
 plot(yt,zy);
 plot(mt,zm,'k');
 plot(yt,fy,'r');
 
%  plot(yt,zdy);
%  plot(mt,zd3m);
 legend('daily','annual','monthly','trend');
  ylabel('Fremantle Sealevel (mm)')
  set(gca,'xlim',[yt(1) mt(end)]);
  datetick('x','yyyy','keeplimits')
  set(gca,'fontsize',14)
  grid on
  box on



%% save it

[myy mmm mdd]=datevec(mt);

rownames={'year','month','day','monthmean','3month moving avg','trend','detrended monthmean' ,'detrended monthmean 3mo moving avg'};

month_table=table(myy',mmm',mdd',zm',z3m',fm',zdm',zd3m','VariableNames',rownames);
writetable(month_table,'Fremantle_sealevel_month_table.csv');


[yyy ymm ydd]=datevec(yt);
yrownames={'year','month','day','annualmean','trend','detrended annualmean' };
year_table=table(yyy',ymm',ydd',zy',fy',zdy','VariableNames',yrownames);
writetable(year_table,'Fremantle_sealevel_annualmean_table.csv');

