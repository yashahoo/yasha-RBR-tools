function [f, P, Pstat]=autospec(x,srate,dof,tcut,icut,scut,wcut);
% AUTOSPEC.M
% CALCULATES THE ONE-SIDED AUTO-SPECTRA FOLLOWING BIN-AVERAGING 
% METHOD AND DETERMINES THE PARTITIONING OF SPECTRAL ENERGY ACROSS
% THE FREQUENCY BANDS
%
% THE SPECTRUM IS NORMALISED SO THAT: 
%    Variance == sum(x.^2)/n = sum(AutoSpec(:,2))*(m*samp_freq/n)  
%
% THE PRODUCED OUTPUT CONSISTS OF:
%	1) f = frequency axis (Hz)
%	2) P = spectral density (m^2/Hz)
%	3) Pstat = 
%	  [Tpeak Tot_Var Trend_Var Inf_Var Swell_Var Wind_Var Noise_Var
%	    Trend_Per Inf_Per Swell_Per Wind_Per Noise_Per 
%	      Hrms Trend_ht Inf_ht Swell_ht Wind_ht Noise_ht]
clf
% SPECIFY TAPER WINDOW (e.g. BOXCAR, HANN, HAMMING)
n=length(x);           	% Number of data points in time series
w=hanning(n);             	% Hann window
Cf=n/sum(w.^2);             	% Correction factor for tapering
xw=w.*detrend(x);          	% Detrend and taper the time series


% ONE-SIDED PERIODOGRAM
Periodo=periodo(xw,srate);  % Calculate one-sided periodogram

% PARTITIONING OF SPECTRAL ENERGY
tcut=floor(tcut*n/srate);		
icut=floor(icut*n/srate);
scut=floor(scut*n/srate);
wcut=floor(wcut*n/srate);

Norm=Cf*srate/n;			% normalising factor

% CALCULATE VARIANCE
Tot_Var=sum(Periodo(:,2))*Norm;		
Trend_Var=sum(Periodo(1:tcut,2))*Norm;
Inf_Var=sum(Periodo(tcut+1:icut,2))*Norm;
Swell_Var=sum(Periodo(icut+1:scut,2))*Norm;
Wind_Var=sum(Periodo(scut+1:wcut,2))*Norm;
Noise_Var=sum(Periodo(wcut+1:length(Periodo),2))*Norm;

% CALCULATE OF ENERGY
Trend_Per=Trend_Var*100/Tot_Var;
Inf_Per=Inf_Var*100/Tot_Var;
Swell_Per=Swell_Var*100/Tot_Var;
Wind_Per=Wind_Var*100/Tot_Var;
Noise_Per=Noise_Var*100/Tot_Var;

% CALCULATE rms WAVE HEIGHT
Hrms=sqrt(8*Tot_Var);	
Trend_ht=sqrt(8*Trend_Var);
Inf_ht=sqrt(8*Inf_Var);
Swell_ht=sqrt(8*Swell_Var);
Wind_ht=sqrt(8*Wind_Var);
Noise_ht=sqrt(8*Noise_Var);

% ONE-SIDED DENSITY SPECTRUM
numbin=dof/2;			% number of estimates in bin
% numest=n/2/numbin;		% Number of spectral estimates
numest=floor(n/2/numbin); % yh try to fix for failure when not using 2^ sample points
P=zeros(numest,2);		% Pre-allocate matrix for spectrum

  for k=0:numest-1           	% Bin counter
    P(k+1,1)=mean(Periodo((k*numbin)+1:(k+1)*numbin,1));
				% Average frequency across bins
    P(k+1,2)=mean(Periodo((k*numbin)+1:(k+1)*numbin,2));
				% Average density across bins
  end
  
P(:,2)=P(:,2)*Cf;		% correct for tapering

% PRODUCE OUTPUT
f=P(:,1);			% frequency
P=P(:,2);			% spectral density
[a,b]=max(P);
Tpeak=1/f(b);			% spectral peak period

Pstat=[Tpeak Tot_Var Trend_Var Inf_Var Swell_Var Wind_Var Noise_Var ...
        Trend_Per Inf_Per Swell_Per Wind_Per Noise_Per ...
          Hrms Trend_ht Inf_ht Swell_ht Wind_ht Noise_ht];
Pstat=Pstat';
