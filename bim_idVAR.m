%% General linear regression model
% given data matrix and indexes of predicted and predictor variables,
% builds observation matrices and performs linear regression through least squares model identification

%%% INPUT
% Y: original time series (N observations x 2 processes) 
% jv: indexes of predicted series
% iv: indexes of predictors
% iv_lags: lags of predictors, in a row vector (they are all the same for the various predictors)

function out=bim_idVAR(Y,jv,iv,iv_lags)

N=size(Y,1);

% compute the maxlag
maxlag=max(iv_lags);

My=Y((maxlag+1:N)',jv); % observation matrix of the predicted variables
MX=[]; % observation matrix of the predictors
for l=1:length(iv_lags)
    for k=1:length(iv)
        MX=[MX Y(maxlag+1-iv_lags(l):N-iv_lags(l),iv(k))]; %#ok
    end
end

eA=(MX'*MX)\(MX'*My); % model coefficients

eu=My-MX*eA; % residuals
es2u=cov(eu); % covariance of residuals
es2y=cov(My); % covariance of predicted variables
erho2xy=1-det(es2u)/det(es2y); % squared correlation

%%% OUTPUT
out.eA=eA; 
out.eu=eu;
out.es2u=es2u;
out.es2y=es2y;
out.erho2=erho2xy;
out.My=My;
out.MX=MX;

end
