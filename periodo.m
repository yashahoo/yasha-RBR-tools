function P = periodo(x, srate) 
% P=periodo(x, srate) 
% Calculate one-sided periodogram
%  INPUT: one column time series (x)
%	      sample frequency in hertz (srate)
%
% The periodogram (P) is normalized so that (length of time series, n):
%  variance of the (detrended) time series = sum of the spectral estimates.
%  The binwidth of the last (Nyquist) frequency bin is half that of the others.
%
% ********************Parseval's Theorum******************************
%    Variance == sum(x.^2)/n = sum(P(:,2))*(srate/n)  
%
%  OUTPUT  P=[Freq Density Phase]
%    Freq       - Hz (midpoint of frequency bin}
%    Density    - units^2/Hz
%    Phase      - radians

%       Authors: Bruce Hegge / Gerd Masselink (17/01/95)
%       Department of Geography / Center for Water Research
%       University of Western Australia
%       Nedlands,  6009
%       bruce@gis.uwa.edu.au / masselin@cwr.uwa.edu.au

% TWO-SIDED PERIODOGRAM
n=length(x);                    	% Number of data points in time series
Xx = zeros(n,1);                	% Pre-allocate output vector
Xx = abs(fft(x)).^2;     

% Power estimates of two-sided periodogram
PhaXx = atan(imag(fft(x))./...  	% Phase estimates of two-sided periodogram
	real(fft(x)));

% ONE-SIDED PERIODOGRAM
% Define first half of periodogram
select = [1;ones((n/2)-1,1).*2;1]; 	% Define first half of periodogram
Xx = Xx(1:(n/2)+1).*select;        	% Power estimates of one-sided periodogram
Xx = Xx/(srate*n); 	            	% Normalizing for Parseval's theorum
freq=(((0:n/2)'*srate)/n); 	    	% Frequency axis of one-sided periodogram
P=[freq Xx PhaXx(1:(n/2)+1)];  		% One-sided periodogram 
					% P=[freq, density, phase]
P(1,:)=[];                         	% Remove the DC value

