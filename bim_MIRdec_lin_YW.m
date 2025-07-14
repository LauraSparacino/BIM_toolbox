%% computation of MIR and of its time-domain decomposition and frequency-domain expansion for two scalar processes
%% estimation through resolution of the Yule Walker equations 

%%% INPUT
% Am=[A(1)...A(p)]: 2*2p matrix of the VAR model coefficients (strictly causal model)
% Su: 2*2 covariance matrix of the input noises
% fs: sampling frequency
% nfft: number of points of the frequency axis
% q: number of lags used to represent the past states of the processes (only for time domain analysis)
% do_td: 'n' or 'y'

function out = bim_MIRdec_lin_YW(Am,Su,fs,nfft,q,do_td)

narginchk(3,6)
if nargin < 6, do_td='n'; end % default: compute time domain measures from spectral measures
if nargin < 5, q=20; end
if nargin < 4, nfft=1000; end
if nargin < 3, fs=1; end 

% compute the parametric spectral density
[S,H,f] = bim_VARspectra(Am,Su,fs,nfft);

% frequency domain measures (coupling and causality measures)
[f12,f1_2,f2_1,f1o2] = bim_fGC_lin(S,H,Su,nfft);

if do_td=='y'
    % time domain measures
    % obtained from resolution of the YW equations
    [~,SigmaS_S] = bim_LinReg(Am,Su,q,[1 2],[1 2]); % bivariate model, [Y1present Y2present] given [Y1past Y2past]
    [~,Sigma1_12] = bim_LinReg(Am,Su,q,1,[1 2]); % full model, Y1present given [Y1past Y2past]
    [~,Sigma2_12] = bim_LinReg(Am,Su,q,2,[1 2]); % full model, Y2present given [Y1past Y2past]

    [~,Sigma1_1] = bim_LinReg(Am,Su,q,1,1); % restricted model, Y1present given Y1past
    [~,Sigma2_2] = bim_LinReg(Am,Su,q,2,2); % restricted model, Y2present given Y2past

    F2_1=log(Sigma1_1/Sigma1_12); % GC from 2 to 1
    F1_2=log(Sigma2_2/Sigma2_12); % GC from 1 to 2
    F1o2=log((Sigma1_12*Sigma2_12)/det(SigmaS_S)); % IC
    % F12=log((Sigma1_1*Sigma2_2)/det(SigmaS_S)); % TD
    F12=F2_1+F1_2+F1o2; % TD 
else
    % time domain measures
    % obtained as the integrals of the corresponding spectral functions
    F2_1=sum(f2_1)/nfft; % GC from 2 to 1
    F1_2=sum(f1_2)/nfft; % GC from 1 to 2
    F1o2=sum(f1o2)/nfft; % IC
    F12=sum(f12)/nfft; % TD
end

%%% OUTPUT
% information-theoretic measures: I=F/2
out.I12=F12/2; % MIR
out.T1_2=F1_2/2; % TE from 1 to 2
out.T2_1=F2_1/2; % TE from 2 to 1
out.I1o2=F1o2/2; % IT
% spectral measures of coupling and causality
out.f12=f12;
out.f1_2=f1_2;
out.f2_1=f2_1;
out.f1o2=f1o2;
out.S=S; % spectrum
out.f=f; % frequency axis

end
