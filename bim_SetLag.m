%% Sets the vector of indexes for series and lags to be used for conditioning

%%% INPUT
% p: vector of embedding dimensions (one for each series; if 0, the series is excluded)
% tau: vector of embedding delays (one for each series)
% u: vector of propagation times (one for each series)
% zerolag: for each series, 1 if zerolag effect is wanted, 0 if not

%%% OUTPUT
% V: list of candidates, dimension Nc*2, Nc is number of candidates; 1st column: index of the signal; 2nd column: index of the lag

function V=bim_SetLag(p,tau,u,zerolag)

M=length(p);
if nargin<4, zerolag=zeros(1,M); end
if nargin<3, u=ones(1,M); end

% Set time series and lags
V=[];
for m=1:M        
    if zerolag(m)==1
        V=[V; [m 0]]; %#ok
    end
    for k=1:p(m)
        V=[V; [m u(m)+tau(m)*(k-1)]];  %#ok
    end
end
