function [PSD] = oppsd(t,serie,window,p1,p2,p3,CL,PPLOT)
% out=oppsd(sampleperiod,data,1,1,64,1.095,1,1);  
% POWER SPECTRAL DENSITY 
% ______________________________________________________________
%
% COMMAND: OpPSD(t,serie,window,p1,p2,p3,CL);
%
% INPUT
% - t:		sampling period;
% - serie: 	time series;
% - window: 	0 no window
%		1 cosine taper window
%		2 Hanning window
% - IF p1 > 0 ==> average in frequency domain (minB,maxB,ratio);
% 		minB = p1: minimum band;
%		maxB = p2: maximum band;
%		ratio= p3: rate of band increase;
% - IF p1 <=0 ==> average in time domain (WSize,WShift);
%		WSize = p2: window size in seconds;
%		WShift= p3: window shift in seconds;
% - CL:		confidency level option (CL=1 means conf. level stored;
% - PPLOT:	0 means no plot & 1 generates plot.
%
%
% OUTPUT:	[Freq PSD]
%
% REFERENCE:
% Bendat, J.S. & Piersol, A.G., 1986. 'Random data: analysis and
%     measurement procedures' 2nd Ed., John Wiley & Sons, 566 pp.
% _____________________________________________________________________

serie = serie(:);
SSize = length(serie);

if (p1 > 0)
  minB = p1;
  maxB = p2;
  ratio = p3;
  w = SSize;
  s = 0;
  TAvg = 1;
  if (PPLOT == 1)
  fprintf('\n  POWER SPECTRAL DENSITY AVERAGED IN FREQUENCY DOMAIN\n\n');
  end;
else
  minB = 1;
  maxB = 1;
  ratio = 1;
  WSize = p2;
  WShift = p3;
  w = WSize/t;
  s = WShift/t;
  TAvg = fix((SSize-w+s)/s);
  if (PPLOT == 1)
  fprintf('\n  POWER SPECTRAL DENSITY AVERAGED IN TIME DOMAIN\n\n');
  end;
end;

index = [1:w];
PSD = zeros(fix(w/2),1);
for i = 1:TAvg
  dserie = (serie(index) - mean(serie(index)));
  %dserie = detrend(serie(index));
  index = index + s;
  if (window == 2)
    var0 = std(dserie)^2;
    h = .5*(1-cos(2*pi*(1:w)'/(w+1)));
    dserie = h.*dserie;
    var1 = std(dserie)^2;
    factor = var0/var1;
    if (PPLOT == 1)
    fprintf('   HANNING WINDOW: variance correction factor is %5.2f \n',factor);
    end;
  elseif (window == 1)
    var0 = std(dserie)^2;
    h = ones(w,1);
    alfa = ceil(0.10*w);
    h(1:alfa) = 0.5*(1 - cos(pi*(1:alfa)/alfa));
    h((w-alfa):w) = 0.5*(1 - cos(pi*((w-alfa):w)/alfa));
    dserie = h.*dserie;
    var1 = std(dserie)^2;
    factor = var0/var1;
    if (PPLOT == 1)
    fprintf('   COSINE TAPER WINDOW: variance correction factor is 5.2f\n',factor);
    end;
else
    factor = 1;
end;

X = fft(dserie);
X = X(2:1+(w/2));
PSD = PSD + (factor*(abs(X).^2)/(w/t));

end;

PSD = PSD/TAvg;
Freq = [1:(w/2)]'/(w*t);

size(Freq);
size(PSD);
if (p1 > 0) & (maxB > 1.0)
  PSD = opsmooth(PSD,minB,maxB,ratio);
  PSD = PSD(:,1);
  Freq = opsmooth(Freq,minB,maxB,ratio);
  Band = Freq(:,2);
  Freq = Freq(:,1);
else
  Band = TAvg*ones(size(Freq));
end;


if (CL == -1) & (PPLOT ==1)
  PSD = opclevel([PSD Band]);
  plot(Freq,PSD(:,1),'k',Freq,PSD(:,2),':',Freq,PSD(:,3),':');
  axis([0 5 0.0001 100]);
  set(gca,'yscale','log');
  xlabel('Frequency');
  ylabel('Spectral Density');

elseif (CL == +1) & (PPLOT == 1)
  CLV = min(PSD)*ones(size(PSD));
  CLV = opclevel([ CLV Band]);
  plot(Freq,PSD(:,1),'k',Freq,CLV(:,2),':r',Freq,CLV(:,3),':r');
  axis([0 5 0.0001 100]);
  set(gca,'yscale','log');
  xlabel('Frequency');
  ylabel('Spectral Density');
else 
  if (PPLOT == 1)
    plot(Freq,PSD(:,1),'k');
    axis([0 5 0.0001 100]);
    set(gca,'yscale','log');
    xlabel('Frequency');
    ylabel('Spectral Density');
  end;
end;

PSD=[Freq PSD(:,1)];
return




