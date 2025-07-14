%% Model order selection for identification of a strictly causal VAR model

%%% INPUT
% Y, 2*N matrix of time series (each time series is in a row) - it works with M > 2 processes too
% pmax, maximum tested model order

%%% OUTPUT
% pottaic: model order optimized with multichannel Akaike Information Criterion (AIC)
% pottaic: model order optimized with  Minimum Description Length (MDL) criterion
% aic: values of AIC index as a function of the model order
% mdl: values of MDL index as a function of the model order

function [pottaic,pottmdl,aic,mdl] = bim_mos_idVAR(Y,pmax)

N=size(Y,2);
M=size(Y,1); 

% figures of merit
aic=NaN*ones(pmax,1); mdl=aic;

for p=1:pmax
    
    out=bim_idVAR(Y',1:M,1:M,1:p);
    Su=out.es2u; % covariance matrix of the residuals

    % multivariate AIC 
    aic(p)=N*log(det(Su))+2*M*M*p; 
    
    % multivariate MDL
    mdl(p)=N*log(det(Su))+log(N)*M*M*p; 
        
end

pottaic=find(aic == min(aic));
pottmdl=find(mdl == min(mdl));
