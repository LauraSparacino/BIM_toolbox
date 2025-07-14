%% computation of MIR and of its time-domain decomposition and frequency-domain expansion for two scalar processes
%% estimation through state space models

%%% INPUT
% Am=[A(1)...A(p)]: 2*2p matrix of the VAR model coefficients (strictly causal model)
% Su: 2*2 covariance matrix of the input noises
% fs: sampling frequency
% nfft: number of points of the frequency axis
% do_td: 'n' or 'y'

function out=bim_MIRdec_lin_SS(Am,Su,fs,nfft,do_td)

narginchk(2,5)
if nargin < 5, do_td='n'; end % default: compute time domain measures from spectral measures
if nargin < 4, nfft=1000; end
if nargin < 3, fs=1; end 

% VAR spectra (from full model)
[S,H,f] = bim_VARspectra(Am,Su,fs,nfft);

% frequency domain measures
[f12,f1_2,f2_1,f1o2] = bim_fGC_lin(S,H,Su,nfft);

if do_td=='y'
    
    % state space model
    [A,C,K]=bim_SSmodel(Am);
    
    % time domain measures
    % obtained from identification of the SS reduced model corresponding to the VAR model
    [~,~,VR1]=bim_submodel(A,C,K,Su,1); % reduced model with the 1st process only 
    F2_1 = log(det(VR1)) - log(det(Su(1,1))); % GC from 2 to 1

    [~,~,VR2]=bim_submodel(A,C,K,Su,2); % reduced model with the 2nd process only 
    F1_2 = log(det(VR2)) - log(det(Su(2,2))); % GC from 1 to 2

    F1o2=log((Su(1,1)*Su(2,2))/det(Su)); % IC
    % F12=log((VR1*VR2)/det(Su)); % TD
    F12=F1_2+F2_1+F1o2; % TD
    
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


