function [bad] = get_outliers_yh(x,nd,wind)
% [badi] = get_outliers_yh(x,nd,wind)
% [input]
% x = data vector
% nd = number of std dev away from median to define threshold
% wind = window definition based on indexes behind/infront of
% centre.
% yasha hetzel 20210518
% based on on this:
% https://www.mathworks.com/help/signal/ug/hampel-filter-algorithm.html
% Use moving window to get indexes of outliers > nd standard deviations
% from the window median. Window  can be balanced or offset around center.
% % eg
% k = [1 3]; % window definition [behind infront]
% nd = 2; % standard deviations threshold

% lx = 24;
% x = randn(1,lx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lx =length(x);
k=wind;
if length(k)==1 % if only single number is given for window it is balanced.
    k=repmat(k,1,2);
    disp(['balanced window used with half size = ' num2str(wind)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define indices of window
iLo = (1:lx)-k(1);
iHi = (1:lx)+k(2);

% Truncate the window so that the function computes medians of smaller segments as it reaches the signal edges.
iLo(iLo<1) = 1;
iHi(iHi>lx) = lx;
% Record the median of each surrounding window. Find the median of the absolute deviation of each element with respect to the window median.
mmed=nan(size(x));mmad=nan(size(x));
for j = 1:length(x)
    w = x(iLo(j):iHi(j));
    medj = median(w,'omitnan');
    mmed(j) = medj;
    mmad(j) = median(abs(w-medj),'omitnan');
end
% Scale the median absolute deviation to obtain an estimate of the standard deviation of a normal distribution.

sd = mmad/(erfinv(1/2)*sqrt(2));
% Find the samples that differ from the median by more than nd = 2 standard deviations. 

bad = abs(x-mmed) > nd*sd;
% badi=find(bad>0);

% Replace each of those outliers by the value of the median of its surrounding window. This is the essence of the Hampel algorithm.

% yu = x;
% yu(ki) = mmed(ki);


end
%% fig
% 
% xx=1:length(x);
% 
% figure;
% plot(x,'k-')
% hold on;
% plot(xx(badi),x(badi),'r.')