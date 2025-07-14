# BIM_toolbox
TOOLBOX DESCRIPTION

LINEAR MODEL-BASED ESTIMATION
•	bim_mos_idVAR: model order selection for identification of a strictly causal vector autoregressive (VAR) model
•	bim_idVAR: general linear regression modelling through least squares model identification
•	bim_VARspectra: VAR spectral matrices (transfer function and power spectral density - PSD)
•	bim_fGC_lin: frequency domain bivariate Granger Causality (GC), Total Dependence (TD) and Instantaneous Causality (IC) from the PSD of a VAR model
Resolution of the Yule-Walker (YW) equations
•	bim_Yule: solution of the YW equations for a VAR process (using discrete time Lyapunov equation)
•	bim_LinReg: linear regression of random processes through resolution of the YW equations; performs linear regression of the present state of given target processes from the past states of given driver processes
•	bim_MIRdec_lin_YW: performs computation of the mutual information rate (MIR) and the causal terms of its decomposition (transfer entropies – TE – and instantaneous transfer – IT) in the time and spectral domains; estimation through resolution of the YW equations
State-space (SS) models
•	bim_SSmodel: computation of SS model parameters [A,C,K] from VAR model parameters [Am]
•	bim_submodel: derivation of a submodel (i.e., a reduced model) of a state space (SS) model 
•	bim_MIRdec_lin_SS: performs computation of the mutual information rate (MIR) and the causal terms of its decomposition (TE, IT) in the time and spectral domains; estimation through SS models

MODEL-FREE ESTIMATION
•	bim_MIRdec_knn: decomposition of the MIR into TEs and IT through the k-nearest neighbors (KNN) estimator
•	bim_MIRdec_bin: decomposition of the MIR into TEs and IT through the binning estimator
•	bim_MIRdec_perm: decomposition of the MIR into TEs and IT through the permutation estimator
•	bim_H: entropy of a discrete multidimensional variable (logarithm of the probability distribution)
•	bim_quantization: quantization of the input series with a given number of quantization levels – used for binning estimator
•	bim_ObsMat: computation of the observation matrix (for entropy computation) 
•	bim_SetLag: sets the vector of indexes for series and lags to be used for conditioning 

OTHER FUNCTIONS
•	bim_AR_filter: autoregressive low-pass, high-pass filter for pre-processing
•	bim_surrtimeshift: time shifted surrogates (makes use of a circular shift with a minimum number of shifted samples)
•	bim_WCspectra: non-parametric power spectral density via weighted covariance estimator

