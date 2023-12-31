function [v,gaps_filled,gaps_unfilled,bind,gapindxj0,gapindxj1]=fillgapYLH(u,maxgap);
% FILLGAP fill gaps (nans) shorter than MAXGAP by linear interpolation
% gaps bigger than MAXGAP values long will not be filled
% Usage:
%  [v,gaps_filled,gaps_unfilled,bind,gapindxj0,gapindxj1]=fillgapYLH(u,maxgap);
% [v,gaps_filled,gaps_unfilled,bind]=fillgap(u,maxgap);
% bind = boundary between good data and bad data with gaps > maxgap
% gapindxj0 / j1 = indexes (start/ end) of gaps > maxgap
% u = vector or 2D matrix (columns of time series data)
% maxgap = max number of values in gap to fill
%
% Note: must be uniformly spaced time series.  If not, use
% "interp1gap" to fill in mixing time and NaN gap values before
% using "fillgap".
% rsignell@usgs.gov

[m,n]=size(u);
u=u(:);
good=find(isfinite(u));
if ~isempty(good)
    dind=diff(good);
    ibound1=find(dind<=(maxgap+1)& dind>1);
    ibound=find(dind>(maxgap+1));
    gaps_filled=length(ibound1);
    gaps_unfilled=length(ibound);
    
    bind=good(ibound); %boundary between good data and bad data with gaps > maxgap
    
    x=1:length(u);
    v=interp_r(x(good),u(good),x);
    v=reshape(v,m,n);
    %
    % mask gaps that were longer than maxgap
    %
    if ~isempty(bind)
        for j=1:length(bind),
            j0=bind(j)+1;
            j1=bind(j)+dind(ibound(j))-1;
            jj=j0:j1;
            v(jj)=v(jj)*NaN;
            gapindxj0(j)=j0;
            gapindxj1(j)=j1;
        end
    else
        gapindxj0=nan;
        gapindxj1 =nan;
    end
    
elseif isempty(good)
    disp(['No data'])
    v=u;
    gaps_filled=nan;
    gaps_unfilled=nan;
    bind=nan;
    gapindxj0=nan;
    gapindxj1 =nan;
    
end

