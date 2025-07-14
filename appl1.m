%% applies MIR decomposition on exemplary RR-SAP time series
%%% linear approach in the time and frequency domains
clear; close all; clc;
addpath([pwd '\..\functions\'])

%%%% parameters

%%% binning estimator 
m_bin=1; % number of past lags of Markov processes
b=4; % n. of bins

%%% permutation estimator 
m_perm=3; % number of past lags of Markov processes

%%% linear estimator
sel_crit='aic'; % 'fix', 'aic', 'mdl'
m_lin_fix=3; % fixed model order
m_lin_max=8; % maximum scanned model order
q=20; % number of lags for estimation of correlations

%%% knn estimator
m_knn=3; 
k=10;

%%% other par
base=0; % 2 for entropy in bits, 0 for entropy in nats
pfilter = 0.94; % high pass filtering
nfft=1000; % number of points of frequency axis
ns=100; % number of surrogates
alpha=0.05; % confidence level
minshift=20; % for time shifted surrogates
range1=[0.04 0.15]; range2=[0.15 0.4]; % frequency bands of interest

%% open data
load([pwd '\data\data1.mat']); % RR[sec]->column1; SAP[mmHg]->column2

Y=data_RR_SAP; % original series

M=size(Y,2); % number of series

mean_RR=mean(Y(:,1)); % [sec]
fs=1/mean_RR; % sampling frequency

% compute bands in frequency bins
nrange1=round((nfft*2/fs)*range1);
if range2(2) < fs/2
    nrange2=round((nfft*2/fs)*range2)+[1 0];
else
    nrange2(1)=round((nfft*2/fs)*range2(1))+1;
    nrange2(2)=round((nfft*2/fs)*fs/2);
end

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

%         Ys(:,1)=bim_surriaafft(Y(:,1));
%         Ys(:,2)=Y(:,2);
        
        Ys(:,1)=bim_surrtimeshift(Y(:,1),minshift); % time shifted surrogates
        Ys(:,2)=Y(:,2);

        data=Ys;
    end
    
    %%% linear approach ( the two methods are equivalent)
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

    % YULE WALKER
    out_lin_YW = bim_MIRdec_lin_YW(Am,Su,fs,nfft,q,'y');
    f=out_lin_YW.f;
    % time domain measures (row 1)
    YW_I12(1,is)=out_lin_YW.I12; % MIR
    YW_T1_2(1,is)=out_lin_YW.T1_2; % TE 1-->2
    YW_T2_1(1,is)=out_lin_YW.T2_1; % TE 2-->1
    YW_I1o2(1,is)=out_lin_YW.I1o2; % IC
    % spectral measures
    i12(:,is)=out_lin_YW.f12/2;
    t1_2(:,is)=out_lin_YW.f1_2/2;
    t2_1(:,is)=out_lin_YW.f2_1/2;
    i1o2(:,is)=out_lin_YW.f1o2/2;
    % integration alongside band B1 (row 2)
    YW_I12(2,is)=sum(i12(nrange1(1):nrange1(2),is))/nfft;
    YW_T1_2(2,is)=sum(t1_2(nrange1(1):nrange1(2),is))/nfft;
    YW_T2_1(2,is)=sum(t2_1(nrange1(1):nrange1(2),is))/nfft;
    YW_I1o2(2,is)=sum(i1o2(nrange1(1):nrange1(2),is))/nfft;
    % integration alongside band B2 (row 3)
    YW_I12(3,is)=sum(i12(nrange2(1):nrange2(2),is))/nfft;
    YW_T1_2(3,is)=sum(t1_2(nrange2(1):nrange2(2),is))/nfft;
    YW_T2_1(3,is)=sum(t2_1(nrange2(1):nrange2(2),is))/nfft;
    YW_I1o2(3,is)=sum(i1o2(nrange2(1):nrange2(2),is))/nfft;

    % STATE SPACE --- it is the same as the method based on the resolution of the YW equations
    out_lin_SS = bim_MIRdec_lin_SS(Am,Su,fs,nfft,'y');
    SS_I12(is)=out_lin_SS.I12; % MIR
    SS_T1_2(is)=out_lin_SS.T1_2; % TE 1-->2
    SS_T2_1(is)=out_lin_SS.T2_1; % TE 2-->1
    SS_I1o2(is)=out_lin_SS.I1o2; % IC

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
% LINEAR (time and frequency domain)
for irange=1:3
    % time
    YW_I12_th(irange)=prctile(YW_I12(irange,2:end),100-alpha*100); % MIR 
    YW_I1_2_th(irange)=prctile(YW_T1_2(irange,2:end),100-alpha*100); % TE 1-->2
    YW_I2_1_th(irange)=prctile(YW_T2_1(irange,2:end),100-alpha*100); % TE 2-->1
    YW_I1o2_th(irange)=prctile(YW_I1o2(irange,2:end),100-alpha*100); % IT 
end
% frequency
i12_th=prctile(i12(:,2:end),[alpha/2 50 100-alpha/2*100],2); % MIR 
t1_2_th=prctile(t1_2(:,2:end),[alpha/2 50 100-alpha/2*100],2); % TE 1-->2
t2_1_th=prctile(t2_1(:,2:end),[alpha/2 50 100-alpha/2*100],2); % TE 2-->1
i1o2_th=prctile(i1o2(:,2:end),[alpha/2 50 100-alpha/2*100],2); % IT

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
disp(['Lin YW: ', num2str(YW_I12(1,1)),' nats']);
disp(['Lin SS: ', num2str(SS_I12(1)),' nats']);
disp(['Knn: ', num2str(k_I12(1)),' nats']);
disp(['Bin: ', num2str(b_I12(1)),' bits']);
disp(['Perm: ', num2str(p_I12(1)),' bits']);
disp(' ')
disp('Estimated values of TE 1-->2:');
disp(['Lin YW: ', num2str(YW_T1_2(1,1)),' nats']);
disp(['Lin SS: ', num2str(SS_T1_2(1)),' nats']);
disp(['Knn: ', num2str(k_T1_2(1)),' nats']);
disp(['Bin: ', num2str(b_T1_2(1)),' bits']);
disp(['Perm: ', num2str(p_T1_2(1)),' bits']);
disp(' ')
disp('Estimated values of TE 2-->1:');
disp(['Lin YW: ', num2str(YW_T2_1(1,1)),' nats']);
disp(['Lin SS: ', num2str(SS_T2_1(1)),' nats']);
disp(['Knn: ', num2str(k_T2_1(1)),' nats']);
disp(['Bin: ', num2str(b_T2_1(1)),' bits']);
disp(['Perm: ', num2str(p_T2_1(1)),' bits']);
disp(' ')
disp('Estimated values of IC:');
disp(['Lin YW: ', num2str(YW_I1o2(1,1)),' nats']);
disp(['Lin SS: ', num2str(SS_I1o2(1)),' nats']);
disp(['Knn: ', num2str(k_I1o2(1)),' nats']);
disp(['Bin: ', num2str(b_I1o2(1)),' bits']);
disp(['Perm: ', num2str(p_I1o2(1)),' bits']);

%% plot results --- only linear method, time and frequency domain 
DimFont=18;
figure('Color','w','WindowState','maximized')

% 1st row: frequency domain distributions
meas{1}=i12_th; meas{2}=t1_2_th; meas{3}=t2_1_th; meas{4}=i1o2_th;
o_meas{1}=i12; o_meas{2}=t1_2; o_meas{3}=t2_1; o_meas{4}=i1o2;
leg={'f_{Y_1;Y_2}','f_{Y_1 \rightarrow Y_2}','f_{Y_2 \rightarrow Y_1}','f_{Y_1 \cdot Y_2}'};
for imeas=1:length(meas)
    subplot(2,4,imeas) 
    area(f,meas{imeas}(:,3),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8],...
        'HandleVisibility','off'); hold on; % area under the 97.5° prctile
    area(f,meas{imeas}(:,1),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8],...
        'HandleVisibility','off'); hold on; % area under the 2.5° prctile
    plot(f,meas{imeas}(:,1),'Color',[0.2 0.2 0.2],'LineWidth',1.2,...
        'HandleVisibility','off'); hold on; % 2.5° prctile
    plot(f,meas{imeas}(:,2),'Color',[0.5 0.5 0.5],'LineWidth',1.2,...
        'HandleVisibility','off'); hold on; % 50° prctile
    plot(f,meas{imeas}(:,3),'Color',[0.2 0.2 0.2],'LineWidth',1.2,...
        'HandleVisibility','off'); hold on; % 97.5° prctile
    plot(f,o_meas{imeas}(:,1),'k','LineWidth',2); hold on; % original
    xlim([0 fs/2]);
    ylim([-0.5 max(o_meas{1}(:,1))+0.1]);
    xticks(0:0.1:fs/2);
    xlabel('f[Hz]');
    if imeas==1
        ylabel('[nats/Hz]')
    end
    legend(leg{imeas}); legend box off
    pbaspect([1 1 1])
end
clear meas 

% 2nd row: time domain distributions (TIME, B1, B2) - LINEAR
x{1}=1; x{2}=3; x{3}=5;
o_meas{1,1}=YW_I12(1,1); o_meas{2,1}=YW_T1_2(1,1); o_meas{3,1}=YW_T2_1(1,1); o_meas{4,1}=YW_I1o2(1,1); % TIME
o_meas{1,2}=YW_I12(2,1); o_meas{2,2}=YW_T1_2(2,1); o_meas{3,2}=YW_T2_1(2,1); o_meas{4,2}=YW_I1o2(2,1); % LF
o_meas{1,3}=YW_I12(3,1); o_meas{2,3}=YW_T1_2(3,1); o_meas{3,3}=YW_T2_1(3,1); o_meas{4,3}=YW_I1o2(3,1); % HF
meas{1,1}=YW_I12_th(1); meas{2,1}=YW_I1_2_th(1); meas{3,1}=YW_I2_1_th(1); meas{4,1}=YW_I1o2_th(1); % TIME
meas{1,2}=YW_I12_th(2); meas{2,2}=YW_I1_2_th(2); meas{3,2}=YW_I2_1_th(2); meas{4,2}=YW_I1o2_th(2); % LF
meas{1,3}=YW_I12_th(3); meas{2,3}=YW_I1_2_th(3); meas{3,3}=YW_I2_1_th(3); meas{4,3}=YW_I1o2_th(3); % HF
tit={'Mutual Information Rate','Transfer Entropy Y_1 \rightarrow Y_2',...
    'Transfer Entropy Y_2 \rightarrow Y_1','Instantaneous Transfer'};
for imeas=1:size(o_meas,1)
    subplot(2,4,4+imeas)
    for irange=1:3
        if o_meas{imeas,irange} > meas{imeas,irange}
            bar(x{irange},o_meas{imeas,irange},'EdgeColor',[0 0 0],'FaceColor',[0 0 0]); hold on;
        else
            bar(x{irange},o_meas{imeas,irange},'EdgeColor',[0 0 0],'FaceColor',[0.8 0.8 0.8]); hold on;
        end
    end
    xlim([x{1}-1 x{3}+1]);
    ylim([-0.01 max(o_meas{1,1})+0.05])
    xticks([x{1} x{2} x{3}]);
    xticklabels({'LIN_{T}','LIN_{B_1}','LIN_{B_2}'});
    if imeas==1
        ylabel('[nats]')
    end
    title(tit{imeas});
    pbaspect([1 1 1])
end
clear meas

% 3rd row: time domain distributions - KNN, PERM, BIN
% o_meas{1,1}=k_I12(1); o_meas{2,1}=k_T1_2(1); o_meas{3,1}=k_T2_1(1); o_meas{4,1}=k_I1o2(1); % KNN
% o_meas{1,2}=p_I12(1); o_meas{2,2}=p_T1_2(1); o_meas{3,2}=p_T2_1(1); o_meas{4,2}=p_I1o2(1); % PERM
% o_meas{1,3}=b_I12(1); o_meas{2,3}=b_T1_2(1); o_meas{3,3}=b_T2_1(1); o_meas{4,3}=b_I1o2(1); % BIN
% meas{1,1}=k_I12_th; meas{2,1}=k_I1_2_th; meas{3,1}=k_I2_1_th; meas{4,1}=k_I1o2_th; % KNN
% meas{1,2}=p_I12_th; meas{2,2}=p_I1_2_th; meas{3,2}=p_I2_1_th; meas{4,2}=p_I1o2_th; % PERM
% meas{1,3}=b_I12_th; meas{2,3}=b_I1_2_th; meas{3,3}=b_I2_1_th; meas{4,3}=b_I1o2_th; % BIN
% for imeas=1:size(o_meas,1)
%     subplot(3,4,8+imeas)
%     for imethod=1:3
%         if o_meas{imeas,imethod} > meas{imeas,imethod}
%             bar(x{imethod},o_meas{imeas,imethod},'EdgeColor',[0 0 0],'FaceColor',[0 0 0]); hold on;
%         else
%             bar(x{imethod},o_meas{imeas,imethod},'EdgeColor',[0 0 0],'FaceColor',[0.8 0.8 0.8]); hold on;
%         end
%     end
%     xlim([x{1}-1 x{3}+1]);
%     ylim([0 o_meas{1,2}+0.1])
%     xticks([x{1} x{2} x{3}]);
%     xticklabels({'KNN','PERM','BIN'});
%     if imeas==1
%         ylabel('[nats]')
%     end
% end

exportgraphics(gcf,[pwd '\..\figures\appl1_raw.pdf'],...
            'Resolution',600,'ContentType','vector');

