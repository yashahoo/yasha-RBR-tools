function out = newdatestr(T)
% This function converts a time, T, from the DATENUM format into
% yyyy-mm-dd HH:MM:SS.SSSSSS
%
% Here the seconds part retains 6 digits past the decimal.  To change 
% this number, change the number after the decimal after
% the last percent in the format string "f".
[y,mo,d,h,min,s] = datevecmx(T, 2);
f = '%4d-%02d-%02d %2d:%2d:%2.2f\n';
out = sprintf(f, [y, mo, d, h, min, s]');