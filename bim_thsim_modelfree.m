%% Sensivity Analysis of model free estimators --- Simulation of cardiovascular interactions
clear; close all; clc;
addpath([pwd '\..\functions\'])

%% Sensitivity of the KNN estimator to k (number of neighbors)
nsim=100; % number of simulations

%%% knn estimator
m_knn=2; 
k=[5 10 15 20 30];

%%% linear estimator
sel_crit='fix'; % 'fix', 'aic', 'bic'
m_lin_fix=2; % fixed model order
m_lin_max=8; % maximum scanned model order
q=20; % number of lags for estimation of correlations

%%% other par
base=0; % 2 for entropy in bits, 0 for entropy in nats
nfft=1000; % number of points of frequency axis
fs=1;

coup=0:0.1:1; % coupling parameter
M=2; % number of series
tau=ones(1,M); % embedding lag always unitary (for knn, bin & perm approaches)
N = 1000; % time series length

% simulate the dynamics
par.poles{1}=([0.8 0.1]); % Oscillations RR X
par.poles{2}=([0.9 0.3]); % Oscillation RESP y
par.Su=[1 1];

hw1=waitbar(0,'Coupling...');
hw2=waitbar(0,'Number of Neighbors...');
hw3=waitbar(0,'Number of Simulation...');

for ic=1:numel(coup)
    par.coup=[2 1 1 1-coup(ic)]; % simulate causal interaction from y to x
    [Am,Su]=bim_theoreticalVAR(M,par); % VAR parameters
    Am=Am';

    %%% theoretical measures 
    out_th = bim_MIRdec_lin_YW(Am,Su,fs,nfft,q);
    % time domain
    T1_2(ic) = out_th.T1_2;
    T2_1(ic) = out_th.T2_1;
    I1o2(ic) =out_th.I1o2;
    I12(ic) = out_th.I12;
    % spectral measures
    f12{ic}=out_th.f12;
    f1_2{ic}=out_th.f1_2;
    f2_1{ic}=out_th.f2_1;
    f1o2{ic}=out_th.f1o2;
    f=out_th.f; % frequency axis
    
    for n_k=1:numel(k)
        for ns=1:nsim
            %%% Realization
            Un = mvnrnd(zeros(1,M),Su,N);
            Yn = bim_VARfilter(Am,Un');
            Yn=Yn';

            %%% Normalization
            Yn = zscore(Yn);

            %%% knn estimation
            out_knn=bim_MIRdec_knn(Yn,m_knn,tau,k(n_k));
            T1_2_knn(ns,ic,n_k) = out_knn.T1_2;
            T2_1_knn(ns,ic,n_k) = out_knn.T2_1;
            I1o2_knn(ns,ic,n_k) = out_knn.I1o2;
            I12_knn(ns,ic,n_k) = out_knn.I12;
            waitbar(ns/nsim,hw3);
        end
        waitbar(n_k/numel(k),hw2);
    end
    waitbar(ic/numel(coup),hw1);
end

hw1.delete();
hw2.delete();
hw3.delete();

%% Plots
figure('WindowState','maximized');
sgtitle('KNN Estimator','FontWeight','bold');
subplot(2,2,1);
plot(coup,I12,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(k)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(k)/2)) * 0.02;

    % Compute statistics
    y_med = median(I12_knn(:,:,n_k));
    y_std = std(I12_knn(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['k = ' num2str(k(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('I_{X;Y} [nats]');
legend('show');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

subplot(2,2,2);
plot(coup,T1_2,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(k)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(k)/2)) * 0.02;

    % Compute statistics
    y_med = median(T1_2_knn(:,:,n_k));
    y_std = std(T1_2_knn(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['k = ' num2str(k(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('T_{X \rightarrow Y} [nats]');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

subplot(2,2,3);
plot(coup,T2_1,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(k)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(k)/2)) * 0.02;

    % Compute statistics
    y_med = median(T2_1_knn(:,:,n_k));
    y_std = std(T2_1_knn(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['k = ' num2str(k(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('T_{Y \rightarrow X} [nats]');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

subplot(2,2,4);
plot(coup,I1o2,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(k)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(k)/2)) * 0.02;

    % Compute statistics
    y_med = median(I1o2_knn(:,:,n_k));
    y_std = std(I1o2_knn(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['k = ' num2str(k(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('I_{X \cdot Y} [nats]');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

%% Sensitivity of the permutation approach to the embedding dimension q
m_perm=2:5; % number of past lags

hw1=waitbar(0,'Coupling...');
hw2=waitbar(0,'Number of Lags...');
hw3=waitbar(0,'Number of Simulation...');
for ic=1:numel(coup)
    par.coup=[2 1 1 1-coup(ic)];
    [Am,Su]=bim_theoreticalVAR(M,par); % VAR parameters
    Am=Am';
    
    for n_k=1:numel(m_perm)
        for ns=1:nsim
            %%% Realization
            Un = mvnrnd(zeros(1,M),Su,N);
            Yn = bim_VARfilter(Am,Un');
            Yn=Yn';
            
            %%% Normalization
            Yn = zscore(Yn);

            %%% perm estimation
            out_perm=bim_MIRdec_perm(Yn,m_perm(n_k),tau,base);
            T1_2_perm(ns,ic,n_k) = out_perm.T1_2;
            T2_1_perm(ns,ic,n_k) = out_perm.T2_1;
            I1o2_perm(ns,ic,n_k) = out_perm.I1o2;
            I12_perm(ns,ic,n_k) = out_perm.I12;

            waitbar(ns/nsim,hw3);
        end
        waitbar(n_k/numel(m_perm),hw2);
    end
    waitbar(ic/numel(coup),hw1);
end
hw1.delete();
hw2.delete();
hw3.delete();

figure('WindowState','maximized');
sgtitle('Permutation Estimator','FontWeight','bold')
subplot(2,2,1);
plot(coup,I12,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(m_perm)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(m_perm)/2)) * 0.02;

    % Compute statistics
    y_med = median(I12_perm(:,:,n_k));
    y_std = std(I12_perm(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['q = ' num2str(m_perm(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('I_{X;Y} [nats]');
legend('show');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

subplot(2,2,2);
plot(coup,T1_2,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(m_perm)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(m_perm)/2)) * 0.02;

    % Compute statistics
    y_med = median(T1_2_perm(:,:,n_k));
    y_std = std(T1_2_perm(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['q = ' num2str(m_perm(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('T_{X \rigtharrow Y} [nats]');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

subplot(2,2,3);
plot(coup,T2_1,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(m_perm)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(m_perm)/2)) * 0.02;

    % Compute statistics
    y_med = median(T2_1_perm(:,:,n_k));
    y_std = std(T2_1_perm(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['q = ' num2str(m_perm(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('T_{Y \rightarrow X} [nats]');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

subplot(2,2,4);
plot(coup,I1o2,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(m_perm)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(m_perm)/2)) * 0.02;

    % Compute statistics
    y_med = median(I1o2_perm(:,:,n_k));
    y_std = std(I1o2_perm(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['q = ' num2str(m_perm(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('I_{X \cdot Y} [nats]');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

%% Sensitivity of binning approach to the number of bins b
m_bin=2; % number of past lags of Markov processes
b=[2 4 6 8 10]; % n. of bins

hw1=waitbar(0,'Coupling...');
hw2=waitbar(0,'Number of Lags...');
hw3=waitbar(0,'Number of Simulation...');
for ic=1:numel(coup)
    par.coup=[2 1 1 1-coup(ic)];
    [Am,Su]=bim_theoreticalVAR(M,par); % VAR parameters
    Am=Am';
    
    for n_k=1:numel(b)
        for ns=1:nsim
            
            %%% Realization
            Un = mvnrnd(zeros(1,M),Su,N);
            Yn = bim_VARfilter(Am,Un');
            Yn=Yn';
            
            %%% Normalization
            Yn = zscore(Yn);

            out_bin=bim_MIRdec_bin(Yn,b(n_k),m_bin,tau,base);
            T1_2_bin(ns,ic,n_k) = out_bin.T1_2;
            T2_1_bin(ns,ic,n_k) = out_bin.T2_1;
            I1o2_bin(ns,ic,n_k) = out_bin.I1o2;
            I12_bin(ns,ic,n_k) = out_bin.I12;
         

            waitbar(ns/nsim,hw3);
        end
        waitbar(n_k/numel(b),hw2);
    end
    waitbar(ic/numel(coup),hw1);
end
hw1.delete();
hw2.delete();
hw3.delete();

% Plot of the analysis
figure('WindowState','maximized');
sgtitle('Binning Estimator','FontWeight','bold')
subplot(2,2,1);
plot(coup,I12,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(b)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(b)/2)) * 0.02;

    % Compute statistics
    y_med = median(I12_bin(:,:,n_k));
    y_std = std(I12_bin(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['b = ' num2str(b(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('I_{X;Y} [nats]');
legend('show');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

subplot(2,2,2);
plot(coup,T1_2,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(b)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(b)/2)) * 0.02;

    % Compute statistics
    y_med = median(T1_2_bin(:,:,n_k));
    y_std = std(T1_2_bin(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['b = ' num2str(b(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('T_{X \rightarrow Y} [nats]');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

subplot(2,2,3);
plot(coup,T2_1,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(b)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(b)/2)) * 0.02;

    % Compute statistics
    y_med = median(T2_1_bin(:,:,n_k));
    y_std = std(T2_1_bin(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['b = ' num2str(b(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('I_{Y \rightarrow X} [nats]');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

subplot(2,2,4);
plot(coup,I1o2,'Color','k','LineWidth',1.5,'DisplayName','Theoretical Value');
for n_k = 1:numel(b)
    hold on;
    % Compute offset x positions
    x_offset = coup + (n_k - ceil(numel(b)/2)) * 0.02;

    % Compute statistics
    y_med = median(I1o2_bin(:,:,n_k));
    y_std = std(I1o2_bin(:,:,n_k), 0, 1); % use std instead of std2 for vectors

    % Plot with offset
    errorbar(x_offset, y_med, y_std, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['b = ' num2str(b(n_k))]);
end
xlim([-0.1 1.1]);
xlabel('Coupling Strength');
ylabel('I_{X \cdot Y} [nats]');
grid on;
ax=gca;
ax.FontSize=16;
ax.LineWidth=2;

% save Results
save Sensitivity_Analysis.mat
