%% Time shifted surrogates
% it explores the significance of short-term correlations using a minimum
% lon-term time shift; it maintains the individual properties of the single series

function [xs,tauvett]=bim_surrtimeshift(x,taumin,numsurr)

narginchk(1,3);
if nargin < 3, numsurr=1; end % default: 1 surrogate series
if nargin < 2, taumin=1; end % default: shift of one sample

N=length(x);
taumax=N-taumin;

tauvett=randperm(taumax-taumin+1)+taumin-1;

xs=zeros(N,numsurr);
for i=1:numsurr
    xs(:,i)=circshift(x,tauvett(i));
end

