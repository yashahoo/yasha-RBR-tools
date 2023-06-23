function wavepar=wavepar_yh(data,srate,icut);
% WAVEPAR.M
% CALCULATES WAVE PARAMETERS ACCORDING USING ZERO-DOWN-CROSSING APPROACH
% [recommended by the International Association of Hydraulic Research]
% THE TIME SERIES IS HIGHPASS FILTERED PRIOR TO ANALYSIS
% wavepar=wavepar(data,srate,icut);
% srate sampling frequency in Hz
% THE OUTPUTS ARE:
%	c     - number of waves
%   Tz    - mean wave period
%	Hrms  - rms wave height
%	Hmax  - largest crest-to-trough wave
%	Tmax  - period of Hmax
%	Hs    - significant wave height; mean height of largest 1/3 of waves
%	Ts    - the mean period of the largest 1/3 of the waves
%   H10   - mean height of largest 1/10 of waves
%	T10   - the mean period of the largest 1/10 of the waves
	
% DATA PREPARATION
data=data-mean(data,'omitnan');		% remove the mean
data=highpass(data,srate,icut);	% highpass filter

% COUNT THE NUMBER OF ZERO DOWN-CROSSINGS
n=length(data);	
c=0;				% wave counter
for i=1:n-1
  if (data(i)>0) ...
           & (data(i+1)<=0),    % check for zero down-crossing
    c=c+1;			% increment wave counter
    w(c)=i;			% update vector with zero down-crossings
  end
end

% CALCULATE MEAN WAVE PERIOD
totaltime=(w(c)-w(1))/srate;
Tz=totaltime/(c-1);		% mean wave period

% DETERMINE CREST, TROUGH, HEIGHT AND PERIOD OF EACH WAVE
for i=1:c-1			
  wv=data((w(i)+1):w(i+1));	% take segment of the data covering 1 wave
  crest(i)=max(wv);		% determine crest of that wave
  trough(i)=abs(min(wv));	% determine trough of that wave
  
%   %yh save indexes too
%   [crest(i),icrest(i)]=max(wv);	
%   [trough(i),itrough(i)]=min(wv);
%   trough(i)=abs(trough(i));
  
  ht(i)=crest(i)+trough(i);	% determine height of that wave
  per(i)=(w(i+1)-w(i))/srate;	% determine period of that wave
end

% DETERMINE THE ROOT MEAN SQUARE WAVE HEIGHT
ht_sqr=ht.^2;			% create vector with squared wave heights
Hrms=sqrt(mean(ht_sqr));	% determine rms wave height

% DETERMINE HEIGHT AND PERIOD OF THE LARGEST WAVE
[a b]=max(ht);			
Hmax=a;				% wave height of largest wave
Tmax=per(b);			% wave period of largest wave

% DETERMINE WAVE HEIGHT AND PERIOD OF 33% OF THE HIGHEST WAVES
[ht_s,ind] = sort(ht);		% sort the height vector + add index vector
co_3=floor(2*(c-1)/3);		% determine cutoff wave height (33%)
ht_3=ht_s(co_3:c-1);		% select only the highest waves
Hs=mean(ht_3);			% calculate 33% wave height
ind_3=ind(co_3:c-1);		% vector with rank numbers of 33% waves
lenind_3=length(ind_3);		% number of 33% waves
hper_3=0;
for i=1:lenind_3
  hper_3=hper_3+per(ind_3(i));		% sum periods of 33% waves
end
Ts=hper_3/lenind_3;		% calculate wave period of 33% waves

% DETERMINE WAVE HEIGHT AND PERIOD OF 10% OF THE HIGHEST WAVES
[ht_s,ind] = sort(ht);		% sort the height vector + add index vector
co_10=floor(9*(c-1)/10);	% determine cutoff wave height (10%)
ht_10=ht_s(co_10:c-1);		% select only the highest waves
H10=mean(ht_10);		% calculate 10% wave height
ind_10=ind(co_10:c-1);		% vector with rank numbers of 10% waves
lenind_10=length(ind_10);	% number of 10% waves
hper_10=0;
for i=1:lenind_10
  hper_10=hper_10+per(ind_10(i));	% sum periods of 10% waves
end
T10=hper_10/lenind_10;		% calculate wave period of 10% waves

% % PRODUCE OUTPUT VECTOR
% wavepar=[c Tz Hrms Hmax Tmax Hs Ts H10 T10];
% wavepar=wavepar';

wavepar.c=c;
wavepar.Tz=Tz; 
wavepar.Hrms=Hrms; 
wavepar.Hmax=Hmax; 
wavepar.Tmax=Tmax;
wavepar.Hs =Hs;
wavepar.Ts =Ts;
wavepar.H10 =H10;
wavepar.T10=T10;
wavepar.crest=crest ;
% wavepar.icrest=icrest;
wavepar.trough=trough;
% wavepar.itrough=itrough;
wavepar.ht=ht;
wavepar.w=w;
