%% Computation of the mutual information rate (MIR) through the permutation estimator
%% and its decomposition into non-causal terms (entropy rates) and causal terms (transfer entropies)

%%% INPUT
% Y: bivariate time series (N x 2)
% m: memory (number of past samples in the embedding vectors)
% tau: vector of embedding delays (one for each series)
% base: base of the logarithm for entropy computation (2, or not pass argument to measure in bits, 0 to measure in nats)

function out=bim_MIRdec_perm(Y,m,tau,base)

M=size(Y,2);

% CE H(Yn|Yn^m)
V_Y=bim_SetLag([0 m],tau,ones(1,M),[0 0]);
[MyY,~]=bim_ObsMat(Y,2,V_Y);
[~,MyYp]=sort(MyY,2); [~,MyYp]=sort(MyYp,2);
[~,MYp]=sort(MyY(:,2:m+1),2); [~,MYp]=sort(MYp,2);
Myp=MyYp(:,1);
HyY=bim_H([Myp MYp],base);
HY=bim_H(MYp,base);
Hy_Y=HyY-HY;

% CE H(Xn|Xn^m)
V_X=bim_SetLag([m 0],tau,ones(1,M),[0 0]);
[MxX,~]=bim_ObsMat(Y,1,V_X);
[~,MxXp]=sort(MxX,2); [~,MxXp]=sort(MxXp,2); 
[~,MXp]=sort(MxX(:,2:m+1),2); [~,MXp]=sort(MXp,2); 
Mxp=MxXp(:,1);
HxX=bim_H([Mxp MXp],base);
HX=bim_H(MXp,base);
Hx_X=HxX-HX;
 
HyXY=bim_H([Myp MYp MXp],base);
HXY=bim_H([MXp MYp],base);
Hy_XY=HyXY-HXY;

% CE H(Xn|Xn^m,Yn^m) 
HxXY=bim_H([Mxp MXp MYp],base);
Hx_XY=HxXY-HXY;

% Transfer entropy X->Y (1->2)
T1_2=Hy_Y-Hy_XY;

% Transfer entropy Y->X (2->1)
T2_1=Hx_X-Hx_XY;

% instantaneous dependence
HxyXY=bim_H([Mxp Myp MXp MYp],base);
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
