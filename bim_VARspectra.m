%% VAR spectral matrices

%%% INPUT
% Am=[A(1)...A(p)]: 2*2p matrix of the VAR model coefficients (strictly causal model)
% Su: 2*2 covariance matrix of the input noises
% fs: sampling frequency
% nfft: number of points of the frequency axis

%%% OUTPUT
% S: Spectral Matrix 
% H: Tranfer Function Matrix 
% f: frequency vector

function [S,H,f] = bim_VARspectra(Am,Su,fs,nfft)

M = size(Am,1); % number of elements in the system 
p = size(Am,2)/M; % p is the order of the VAR model

if nargin<2, Su = eye(M,M); end % if not specified, we assume uncorrelated noises with unit variance as inputs 
if nargin<3, fs = 1; end
if nargin<4, nfft = 1000; end   
if all(size(nfft)==1)	          % if nfft is scalar
    f = (0:nfft-1)*(fs/(2*nfft)); % frequency axis
else                              % if nfft is a vector, we assume that it is the vector of the frequencies
    f = nfft; nfft = length(nfft);
end

% s = exp(1i*2*pi*f/fs); % vector of complex exponentials
z = 1i*2*pi/fs; 

% Initializations: spectral matrices have M rows, M columns and are calculated at each of the nfft frequencies
H=zeros(M,M,nfft); % Transfer Matrix
S=zeros(M,M,nfft); % Spectral Matrix

A = [eye(M) -Am]; % matrix from which M*M blocks are selected to calculate spectral functions

%% computation of spectral functions
for n=1:nfft % at each frequency
    
    %%% Coefficient matrix in the frequency domain
    As = zeros(M,M); % matrix As(z)=I-sum(A(k))
    for k = 1:p+1
        As = As + A(:,k*M+(1-M:0))*exp(-z*(k-1)*f(n)); % indicization (:,k*M+(1-M:0)) extracts the k-th M*M block from the matrix B (A(1) is in the second block, and so on)
    end

    %%% Transfer matrix 
    H(:,:,n)  = inv(As);

    %%% Spectral matrix 
    S(:,:,n)  = H(:,:,n)*Su*H(:,:,n)'; % ' stands for Hermitian transpose
       
end

