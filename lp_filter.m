function [y,B,A]=lowpass_filter(red,cut,t_step,x)

% funkcija racuna koeficijente za LP butter filter i filtrira niz, 
% poziva se preko argumenata:
% [y,B,A]=function(red,cut,t_step,x) gdje red predstavlja stupanj filtera
% cut je cutoff !!period!!  *u satima* !!!
% t_step je vremenski korak kojim su uzorkovani podaci (recimo 15min) u
% satima!!!! (dakle za dt = 15min, t_step = 0.25 sati) i x je ulazno polje
% koje treba filtrirati filt filt (napred nazad), y je profiltrirano polje
% ivica, 7.6.2003.
f_nyquist=1/t_step/2;
f_cutof=1/cut/f_nyquist;
%disp(['Cut off frekvencija je :'  num2str(f_cutof)]);
[B,A]=butter(red,f_cutof,'low');
y=filtfilt(B,A,x);


