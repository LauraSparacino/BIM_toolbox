%% frequency domain bivariate Granger Causality, Total Dependence and Instantaneous Causality from the PSD of a VAR model

%%% INPUT
% S: spectral matrix - dimension: 2 x 2 x nfft
% H: transfer function matrix - dimension: 2 x 2 x nfft
% Cov: covariance matrix of the residuals of a VAR model (Su) OR covariance matrix of the reduced ISS model (VR) - dimension: 2 x 2
% nfft: number of points for calculation of the spectral functions 

%%% OUTPUT
% f12: Total Dependence (TD)
% f1_2: Granger Causality (GC) from 1 to 2
% f2_1: Granger Causality from 2 to 1
% f1o2: Instantaneous Causality (IC)

function [f12,f1_2,f2_1,f1o2] = bim_fGC_lin(S,H,Cov,nfft)

narginchk(3,4)
if nargin<4, nfft=1000; end

f12=nan*ones(nfft,1); % spectral TD
f1_2=nan*ones(nfft,1); f2_1=f1_2; % Geweke spectral GC from sources to destinations

for n=1:nfft % at each frequency
            
    f2_1(n) = log( abs(S(1,1,n)) / abs(H(1,1,n)*Cov(1,1)*H(1,1,n)'));
    f1_2(n) = log( abs(S(2,2,n)) / abs(H(2,2,n)*Cov(2,2)*H(2,2,n)'));
    
    f12(n)=log( abs(S(1,1,n))*abs(S(2,2,n)) / abs(det(S(:,:,n))) );
                   
end

f1o2=f12-f1_2-f2_1; % spectral IC (from Geweke decomposition)

end


