%% Solution of Yule-Walker Equations for a VAR process (using discrete time Lyapunov equation)

%%% INPUT
% Am  -   generalized connectivity matrix A=(A_1 A_2 ... A_p)
% Su  -   covariance matrix of the residuals
% q   -   number of lags for which to compute correlations

%%% OUTPUT
% R, 2x2x(q+1) matrix of process covariances for lags k=0,1,...,q

function R = bim_Yule(Am,Su,q)

M = size(Am,1); % number of elements in system (M=2)
p=floor(size(Am,2)/M); % number of lags in VAR model

R=NaN*ones(M,M,q+1); % prepare covariance matrices, (:,:,1) is lag 0, (:,:,q+1) is lag q 

% Obtain F and Delta
Im=eye(M*p);
Psi=[Am;Im(1:end-M,:)];
Theta=zeros(p*M,p*M);
Theta(1:M,1:M)=Su(:,:);

% Obtain BigSigma solving the Lyapunov equation: BigSigma = A * BigSigma * A^T + Theta
BigSigma=dlyap(Psi,Theta);

% extract R(0),...,R(p-1)
for i=1:p
    R(:,:,i)=BigSigma(1:M,M*(i-1)+1:M*i);
end

% Yule-Walker solution  for lags >= p
for k=p+1:q+1
    Rk=R(:,:,k-1:-1:k-p);
    Rm=[];
    for ki=1:p
        Rm=[Rm; Rk(:,:,ki)]; %#ok
    end
    R(:,:,k)=Am*Rm;
end



