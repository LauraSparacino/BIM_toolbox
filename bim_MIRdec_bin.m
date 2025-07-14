%% Computation of the mutual information rate (MIR) through the binning estimator
%% and its decomposition into non-causal terms (entropy rates) and causal terms (transfer entropies)

%%% INPUT
% Y: bivariate time series (N x 2)
% b: number of quantization bins
% m: memory (number of past samples in the embedding vectors)
% tau: vector of embedding delays (one for each series)
% base: base of the logarithm for entropy computation (2, or not pass argument to measure in bits, 0 to measure in nats)

function out=bim_MIRdec_bin(Y,b,m,tau,base)

M=size(Y,2);
for im=1:M
    Yq(:,im)=bim_quantization(Y(:,im),b);
end

% conditional entropy H(Yn|Yn^m)
V_Y=bim_SetLag([0 m],tau,ones(1,M),[0 0]);
[MyY,MY]=bim_ObsMat(Yq,2,V_Y);
HyY=bim_H(MyY,base);
HY=bim_H(MY,base);
Hy_Y=HyY-HY;

% conditional entropy H(Yn|Xn^m,Yn^m)
V_XY=bim_SetLag([m m],tau,ones(1,M),[0 0]);
[MyXY,MXY]=bim_ObsMat(Yq,2,V_XY);
HyXY=bim_H(MyXY,base);
HXY=bim_H(MXY,base);
Hy_XY=HyXY-HXY;

% Transfer entropy X->Y (1->2)
T1_2=Hy_Y-Hy_XY;

% conditional entropy H(Xn|Xn^m)
V_X=bim_SetLag([m 0],tau,ones(1,M),[0 0]);
[MxX,MX]=bim_ObsMat(Yq,1,V_X);
HxX=bim_H(MxX,base);
HX=bim_H(MX,base);
Hx_X=HxX-HX;

% conditional entropy H(Xn|Xn^m,Yn^m)
MxXY=bim_ObsMat(Yq,1,V_XY);
HxXY=bim_H(MxXY,base);
Hx_XY=HxXY-HXY;

% Transfer entropy Y->X (2->1)
T2_1=Hx_X-Hx_XY;

% instantaneous dependence
V_XY0=bim_SetLag([m m],tau,ones(1,M),[0 1]);
MxyXY=bim_ObsMat(Yq,1,V_XY0);
HxyXY=bim_H(MxyXY,base);
Hx_yXY=HxyXY-HyXY;
I1o2=Hx_XY-Hx_yXY; 

% Mutual Information Rate
I12=T1_2+T2_1+I1o2;

%%% OUTPUT
% information-theoretic measures
out.I12=I12;  % MIR
out.T1_2=T1_2; % TE 1-->2
out.T2_1=T2_1;  % TE 2-->1
out.I1o2=I1o2;  % IT
