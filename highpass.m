function highpass=highpass(data,srate,cutoff);
% HIGHPASS.M
% HIGH PASS FILTERS A DATA SERIES USING FOURIER TECHNIQUES
% CUT-OFF FREQUENCY IS OPTIONABLE
% srate = sample frequency in Hz



n=length(data);
fftdata=fft(data,n);
cutoff=floor(cutoff/((srate/2)/(n/2)));
fftdata(1:cutoff)=zeros(cutoff,1);
fftdata((n+1)-cutoff:n)=zeros(cutoff,1);
highpass=real(ifft(fftdata));

