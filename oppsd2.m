function [PSD] = oppsd(t,serie,window,p1,p2,p3,CL,PPLOT)

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
  %fprintf('\n  POWER SPECTRAL DENSITY AVERAGED IN FREQUENCY DOMAIN\n\n');
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
    %fprintf('   COSINE TAPER WINDOW: variance correction factor is 5.2f\n',factor);
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
  a1=axes
  set(a1,'xscale','log');
  set(a1,'yscale','log');
  set(a1,'xlim',[1e-07 1e-02]);
  set(a1,'yticklabel',[]);
  set(gca,'fontsize',[14]);
  a2=axes
  loglog(Freq,PSD(:,1),'w',Freq,PSD(:,2),':',Freq,PSD(:,3),':','linewidth',[2]);
  set(gca,'fontsize',[24]);
  ticks=[1e-07,2e-07,3e-07,4e-07,5e-07,6e-07,7e-07,8e-07,9e-07,1e-06,2e-06,3e-06,4e-06,5e-06,6e-06,7e-06,8e-06,9e-06,1e-05,2e-05,3e-05,4e-05,5e-05,6e-05,7e-05,8e-05,9e-05,1e-04,2e-04,3e-04,4e-04,5e-04,6e-04,7e-04,8e-04,9e-04,1e-03,2e-03,3e-03,4e-03,5e-03,6e-03,7e-03,8e-03,9e-03,1e-02,2e-02,3e-02,4e-02,5e-02,6e-02,7e-02,8e-02,9e-02];
  set(a2,'xtick',ticks);
  set(a2,'xticklabel',[]);
  %set(a2,'yticklabel',[]);
  set(a2,'xlim',[1e-07 1e-02]);
  xlabel('Frequency');
  ylabel('Spectral Density');

elseif (CL == +1) & (PPLOT == 1)
  CLV = min(PSD)*ones(size(PSD));
  CLV = opclevel([ CLV Band]);
  a1=axes;
  set(a1,'xscale','log');
  set(a1,'yscale','log');
  set(a1,'xlim',[1e-07 1e-02]);
  set(a1,'yticklabel',[]);
  set(gca,'fontsize',[14]);
  a2=axes;
  loglog(Freq,PSD(:,1),'k',Freq,CLV(:,2),':r',Freq,CLV(:,3),':r','linewidth',[1.5]);
  set(gca,'fontsize',[14]);
  ticks=[1e-07,2e-07,3e-07,4e-07,5e-07,6e-07,7e-07,8e-07,9e-07,1e-06,2e-06,3e-06,4e-06,5e-06,6e-06,7e-06,8e-06,9e-06,1e-05,2e-05,3e-05,4e-05,5e-05,6e-05,7e-05,8e-05,9e-05,1e-04,2e-04,3e-04,4e-04,5e-04,6e-04,7e-04,8e-04,9e-04,1e-03,2e-03,3e-03,4e-03,5e-03,6e-03,7e-03,8e-03,9e-03,1e-02,2e-02,3e-02,4e-02,5e-02,6e-02,7e-02,8e-02,9e-02];
  set(a2,'xtick',ticks);
  set(a2,'xticklabel',[]);
%  set(a2,'yticklabel',[]);
%  set(a2,'xlim',[1e-07 1e-02]);
%  xlabel('Frequency');
%  ylabel('Spectral Density');
else 
  if (PPLOT == 1)
    a1=axes
    set(a1,'xscale','log');
    set(a1,'yscale','log');
    set(a1,'xlim',[1e-07 1e-02]);
    set(a1,'yticklabel',[]);
    set(gca,'fontsize',[20]);
    a2=axes
    loglog(Freq,PSD(:,1),'w','linewidth',[2]);
    set(gca,'fontsize',[24]);
    ticks=[1e-07,2e-07,3e-07,4e-07,5e-07,6e-07,7e-07,8e-07,9e-07,1e-06,2e-06,3e-06,4e-06,5e-06,6e-06,7e-06,8e-06,9e-06,1e-05,2e-05,3e-05,4e-05,5e-05,6e-05,7e-05,8e-05,9e-05,1e-04,2e-04,3e-04,4e-04,5e-04,6e-04,7e-04,8e-04,9e-04,1e-03,2e-03,3e-03,4e-03,5e-03,6e-03,7e-03,8e-03,9e-03,1e-02,2e-02,3e-02,4e-02,5e-02,6e-02,7e-02,8e-02,9e-02];
    set(a2,'xtick',ticks);
    set(a2,'xticklabel',[]);
    %set(a2,'yticklabel',[]);
    set(a2,'xlim',[1e-07 1e-02]);
    xlabel('Frequency');
    ylabel('Spectral Density');
  end;
end;

PSD=[Freq PSD(:,1)];
return
