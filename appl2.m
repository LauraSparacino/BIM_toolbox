%% applies MIR decomposition on exemplary RR-SAP time series
%%% compares binning, permutation, knn and linear approaches
clear; close all; clc;
addpath([pwd '\..\functions\'])

%%%% parameters

%%% binning estimator
m_bin=2; % number of past lags of Markov processes
b=3; % n. of bins

%%% permutation estimator
m_perm=3; % number of past lags of Markov processes

%%% linear estimator
sel_crit='aic'; % 'fix', 'aic', 'mdl'
m_lin_fix=3; % fixed model order
m_lin_max=12; % maximum scanned model order
q=20; % number of lags for estimation of correlations

%%% knn estimator
m_knn=3; 
k=10;

%%% other par
base=0; % 2 for entropy in bits, 0 for entropy in nats
pfilter = 0.94; % high pass filtering
fs=1; % sampling period: 1 month
nfft=1000; % number of points of frequency axis
ns=100; % number of surrogates
alpha=0.05; % confidence level
minshift=20; % time shifted surrogates

%% open data
load([pwd '\data\data2.mat']);
Y=data_clima; % original series 

M=size(Y,2); % number of series

% preprocessing: high-pass filtering
for im = 1:M
    Y(:,im)=bim_AR_filter(Y,im,pfilter);
end

Y=zscore(Y); % normalization (not necessary for bin and perm)

tau=ones(1,M); % embedding lag always unitary (for knn, bin & perm approaches)

%% estimation of bivariate measures of coupling and causality (both original and surrogate data)
for is=1:ns+1
    
    %%% build surrogate data
    if is==1    % original series
        data=Y; 
    else        % realizations of surrogate series
        Ys(:,1)=bim_surrtimeshift(Y(:,1),minshift); % time shifted surrogates
        Ys(:,2)=Y(:,2);

        data=Ys;
    end
    
    %%% linear approach (the two methods are equivalent)
    switch sel_crit
        case 'fix'
            m_lin=m_lin_fix;
        case 'aic'
            m_lin = bim_mos_idVAR(data',m_lin_max); % number of past lags of Markov processes selected through the Akaike Criterion
        case 'mdl'
            [~,m_lin] = bim_mos_idVAR(data',m_lin_max); % number of past lags of Markov processes selected through the Bayesian Criterion
    end
    out=bim_idVAR(data,1:M,1:M,1:m_lin); % VAR model
    Am=out.eA; Am=Am'; % model coefficients
    Su=out.es2u; % covariance matrix of the residuals

    %%% linear (YW)
    out_lin_YW = bim_MIRdec_lin_YW(Am,Su,fs,nfft,q);
    % time domain measures 
    YW_I12(is)=out_lin_YW.I12; % MIR
    YW_T1_2(is)=out_lin_YW.T1_2; % TE 1-->2
    YW_T2_1(is)=out_lin_YW.T2_1; % TE 2-->1
    YW_I1o2(is)=out_lin_YW.I1o2; % IC

    %%% k nearest-neighbor approach
    out_knn = bim_MIRdec_knn(data,m_knn,tau,k);
    k_I12(is)=out_knn.I12;
    k_T1_2(is)=out_knn.T1_2;
    k_T2_1(is)=out_knn.T2_1;
    k_I1o2(is)=out_knn.I1o2;

    %%% binning approach
    out_bin=bim_MIRdec_bin(data,b,m_bin,tau,base);
    b_I12(is)=out_bin.I12;
    b_T1_2(is)=out_bin.T1_2;
    b_T2_1(is)=out_bin.T2_1;
    b_I1o2(is)=out_bin.I1o2;

    %%% permutation approach
    out_perm=bim_MIRdec_perm(data,m_perm,tau,base);
    p_I12(is)=out_perm.I12;
    p_T1_2(is)=out_perm.T1_2;
    p_T2_1(is)=out_perm.T2_1;
    p_I1o2(is)=out_perm.I1o2;

end

%% test significance of coupling and causality measures
% LINEAR 
YW_I12_th=prctile(YW_I12(2:end),100-alpha*100); % MIR 
YW_I1_2_th=prctile(YW_T1_2(2:end),100-alpha*100); % TE 1-->2
YW_I2_1_th=prctile(YW_T2_1(2:end),100-alpha*100); % TE 2-->1
YW_I1o2_th=prctile(YW_I1o2(2:end),100-alpha*100); % IT 

% KNN
k_I12_th=prctile(k_I12(2:end),100-alpha*100); % MIR 
k_I1_2_th=prctile(k_T1_2(2:end),100-alpha*100); % TE 1-->2
k_I2_1_th=prctile(k_T2_1(2:end),100-alpha*100); % TE 2-->1
k_I1o2_th=prctile(k_I1o2(2:end),100-alpha*100); % IT 

% PERMUTATION
p_I12_th=prctile(p_I12(2:end),100-alpha*100); % MIR 
p_I1_2_th=prctile(p_T1_2(2:end),100-alpha*100); % TE 1-->2
p_I2_1_th=prctile(p_T2_1(2:end),100-alpha*100); % TE 2-->1
p_I1o2_th=prctile(p_I1o2(2:end),100-alpha*100); % IT 

% BINNING
b_I12_th=prctile(b_I12(2:end),100-alpha*100); % MIR 
b_I1_2_th=prctile(b_T1_2(2:end),100-alpha*100); % TE 1-->2
b_I2_1_th=prctile(b_T2_1(2:end),100-alpha*100); % TE 2-->1
b_I1o2_th=prctile(b_I1o2(2:end),100-alpha*100); % IT 

%% display (values on the original series)
disp('Estimated values of MIR:');
disp(['Lin YW: ', num2str(YW_I12(1)),' nats']);
disp(['Knn: ', num2str(k_I12(1)),' nats']);
disp(['Bin: ', num2str(b_I12(1)),' bits']);
disp(['Perm: ', num2str(p_I12(1)),' bits']);
disp(' ')
disp('Estimated values of TE 1-->2:');
disp(['Lin YW: ', num2str(YW_T1_2(1)),' nats']);
disp(['Knn: ', num2str(k_T1_2(1)),' nats']);
disp(['Bin: ', num2str(b_T1_2(1)),' bits']);
disp(['Perm: ', num2str(p_T1_2(1)),' bits']);
disp(' ')
disp('Estimated values of TE 2-->1:');
disp(['Lin YW: ', num2str(YW_T2_1(1)),' nats']);
disp(['Knn: ', num2str(k_T2_1(1)),' nats']);
disp(['Bin: ', num2str(b_T2_1(1)),' bits']);
disp(['Perm: ', num2str(p_T2_1(1)),' bits']);
disp(' ')
disp('Estimated values of IC:');
disp(['Lin YW: ', num2str(YW_I1o2(1)),' nats']);
disp(['Knn: ', num2str(k_I1o2(1)),' nats']);
disp(['Bin: ', num2str(b_I1o2(1)),' bits']);
disp(['Perm: ', num2str(p_I1o2(1)),' bits']);

%% plot results
DimFont=18;
figure('Color','w','WindowState','maximized')
tit={'Mutual Information Rate','Transfer Entropy Y_1 \rightarrow Y_2',...
    'Transfer Entropy Y_2 \rightarrow Y_1','Instantaneous Transfer'};

% time domain distributions - LIN, KNN, PERM, BIN
x{1}=1; x{2}=3; x{3}=5; x{4}=7;

o_meas(1,1)=YW_I12(1); o_meas(2,1)=YW_T1_2(1); o_meas(3,1)=YW_T2_1(1); o_meas(4,1)=YW_I1o2(1); % LIN
o_meas(1,2)=k_I12(1); o_meas(2,2)=k_T1_2(1); o_meas(3,2)=k_T2_1(1); o_meas(4,2)=k_I1o2(1); % KNN
o_meas(1,3)=p_I12(1); o_meas(2,3)=p_T1_2(1); o_meas(3,3)=p_T2_1(1); o_meas(4,3)=p_I1o2(1); % PERM
o_meas(1,4)=b_I12(1); o_meas(2,4)=b_T1_2(1); o_meas(3,4)=b_T2_1(1); o_meas(4,4)=b_I1o2(1); % BIN

meas(1,1)=YW_I12_th; meas(2,1)=YW_I1_2_th; meas(3,1)=YW_I2_1_th; meas(4,1)=YW_I1o2_th; % LIN
meas(1,2)=k_I12_th; meas(2,2)=k_I1_2_th; meas(3,2)=k_I2_1_th; meas(4,2)=k_I1o2_th; % KNN
meas(1,3)=p_I12_th; meas(2,3)=p_I1_2_th; meas(3,3)=p_I2_1_th; meas(4,3)=p_I1o2_th; % PERM
meas(1,4)=b_I12_th; meas(2,4)=b_I1_2_th; meas(3,4)=b_I2_1_th; meas(4,4)=b_I1o2_th; % BIN

for imeas=1:size(o_meas,1)
    subplot(2,2,imeas)
    for imethod=1:size(o_meas,2)
        if o_meas(imeas,imethod) > meas(imeas,imethod)
            bar(x{imethod},o_meas(imeas,imethod),'EdgeColor',[0 0 0],'FaceColor',[0 0 0]); hold on;
        else
            bar(x{imethod},o_meas(imeas,imethod),'EdgeColor',[0 0 0],'FaceColor',[0.8 0.8 0.8]); hold on;
        end
    end
    xlim([x{1}-1 x{end}+1]);
    % ylim([0 max(o_meas(imeas,:)) + max(o_meas(imeas,:))/10])
    ylim([0 0.36])
    xticks(x{1}:2:x{end});
    xticklabels({'LIN','KNN','PERM','BIN'});
    ylabel('[nats]')
    title(tit{imeas})
    pbaspect([1 1 1])
end

exportgraphics(gcf,[pwd '\..\figures\appl2_raw.pdf'],...
            'Resolution',600,'ContentType','vector');
        