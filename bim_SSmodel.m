%% computation of state space (SS) model parameters corresponding to a vector autoregressive (VAR) model
% Adapted from SSGC Toolbox: L. Barnett and A. K. Seth, Granger causality for state-space models, Phys.Rev. E 91(4) Rapid Communication, 2015.

%%% INPUT
% Am=[A(1)...A(p)]: 2*2p matrix of the VAR model coefficients

%%% OUTPUT
% A,C,K: innovations form state space parameters

function [A,C,K]=bim_SSmodel(Am)

M=size(Am,1); % number of elements in the system 
p = size(Am,2)/M; % p is VAR model order
pM1 = (p-1)*M;

% A,C,K: innovations form state space parameters
C = Am;
A = [C; eye(pM1) zeros(pM1,M)];
K = [eye(M); zeros(pM1,M)];

% rho: AR spectral norm; rho >= 1 indicates an unstable AR process: rho > 1 is explosive, rho close to 1 may be unit-root
rho = max(abs(eig(A,'nobalance')));
if rho >= 1
    error('unstable AR process')
end

