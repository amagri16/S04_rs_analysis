function [results, dscrp, idx_ft_all] = AD_feature_calculation(y,ch_labels,~,f,filters,filters_PLV,cfg)
f_s = f;
split_len = 500;
filters = filters;
overlap = 0.5;
signal = y;
labels = ch_labels;
results = [];
dscrp = {};
idx_ft_all = [];

%% Feature calculation
%[results, dscrp, idx_ft_all.Pearson] = find_Pearson(signal,labels,250,results,dscrp);
%[results, dscrp, idx_ft_all.Spearman] = find_Spearman(signal,labels,250,results,dscrp);
%[results, dscrp, idx_ft_all.Distcorr] = find_Distcorr(signal,labels,128,results,dscrp);
%[results, dscrp, idx_ft_all.Bandpass_Spearman] = find_bp_Spearman(signal,labels,128,filters,results,dscrp);
%[results, dscrp, idx_ft_all.Bandpass_Distcorr] = find_bp_Distcorr(signal,labels,128,filters,results,dscrp);
%[results, dscrp, idx_ft_all.Cross_Correlation] = find_XCOR(signal,labels,250,results,dscrp);
%[results, dscrp, idx_ft_all.Cross_Covariance] = find_XCOV(signal,labels,250,results,dscrp);
%[results, dscrp, idx_ft_all.PSDs] = find_PSDs(signal,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.FFT] = find_FFT(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.Power_Band] = find_PMTM(signal,f,labels,split_len/2,results,dscrp);
%[results, dscrp, idx_ft_all.Power_Band] = find_PMTM_low_high(signal,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.Coherence] = find_CXY(signal,f,labels,250*6,results,dscrp);
%[results, dscrp, idx_ft_all.DTF_PDC] = find_DTF(signal,f_s,labels,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_Amplitude_Coupling] = find_PAC(signal,filters_PLV,labels,500,results,dscrp);
%[results, dscrp, idx_ft_all.Synchronization_Likelihood] = find_SL(signal,labels,250*2,results,dscrp);
%[results, dscrp, idx_ft_all.Synchronization_Likelihood] = find_SL_bivariate(signal,labels,250*2,results,dscrp);

%[results, dscrp, idx_ft_all.Synchronization_Likelihood] = find_SL_BP(signal,labels,250*2,filters_PLV{1,1},results,dscrp,'delta');
%[results, dscrp, idx_ft_all.Synchronization_Likelihood] = find_SL_BP(signal,labels,250*2,filters_PLV{1,2},results,dscrp,'theta');
%[results, dscrp, idx_ft_all.Synchronization_Likelihood] = find_SL_BP(signal,labels,250*2,filters_PLV{1,3},results,dscrp,'alpha');
%[results, dscrp, idx_ft_all.Synchronization_Likelihood] = find_SL_BP(signal,labels,250*2,filters_PLV{1,4},results,dscrp,'beta');
%[results, dscrp, idx_ft_all.Bispectrum] = find_Bispectrum(signal,f_s,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.Cross_Bispectrum] = find_Cross_Bispectrum(signal,f_s,labels,split_len,results,dscrp);

%[results, dscrp, idx_ft_all.Transfer_Entropy] = find_TE_final(signal,labels,250,results,dscrp);
%[results, dscrp, idx_ft_all.Transfer_Entropy_pink_reduced] = find_TE_pink_reduced(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_Transfer_Entropy] = find_phase_TE(signal,labels,256*6,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_Transfer_Entropy] = find_BP_phase_TE(signal,labels,256*6,results,dscrp,filters_PLV);
%[results, dscrp, idx_ft_all.Phase_Transfer_Entropy] = find_BP_phase_TE_within_ch(signal,labels,256*6,results,dscrp,filters_PLV);
% 
%[results, dscrp, idx_ft_all.Wavelet_Transfer_Entropy] = find_wavelet_TE(signal,f,labels,250,results,dscrp);
%[results, dscrp, idx_ft_all.Wavelet_Transfer_Entropy] = find_wavelet_TE_within_ch(signal,f,labels,250,results,dscrp);
%[results, dscrp, idx_ft_all.BP_Transfer_Entropy] = find_TE_BP(signal,labels,256,results,dscrp,filters);
%[results, dscrp, idx_ft_all.BP_Transfer_Entropy] = find_BP_TE_within_ch(signal,labels,256*6,results,dscrp,filters);

%[results, dscrp, idx_ft_all.Mutual_Information] = find_MI(signal,labels,250,results,dscrp);
%[results, dscrp, idx_ft_all.Temporal_Mutual_Information] = find_TMI(signal,labels,128,results,dscrp);
%[results, dscrp, idx_ft_all.Auto_Mutual_Information] = find_AMI(signal,labels,250,results,dscrp);
 
%[results, dscrp, idx_ft_all.Spectral_Entropy] = find_SE(signal,f_s,labels,split_len,results,dscrp);

%[results, dscrp, idx_ft_all.O_info] = find_O_info(signal,labels,500,6,filters_PLV{1,1},results,dscrp,'delta');
%[results, dscrp, idx_ft_all.O_info] = find_O_info(signal,labels,500,6,filters_PLV{1,2},results,dscrp,'theta');
%[results, dscrp, idx_ft_all.O_info] = find_O_info(signal,labels,500,6,filters_PLV{1,3},results,dscrp,'alpha');
%[results, dscrp, idx_ft_all.O_info] = find_O_info(signal,labels,500,6,filters_PLV{1,4},results,dscrp,'beta');
%[results, dscrp, idx_ft_all.O_info] = find_O_info_channel_permutation(signal,labels,500,6,filters_PLV,results,dscrp,randomPermsWithRep);
% 
%[results, dscrp, idx_ft_all.O_info] = find_O_info(signal,labels,500,6,[],results,dscrp,'');
%[results, dscrp, idx_ft_all.O_info] = find_O_info(signal,labels,500,8,[],results,dscrp,'');
%[results, dscrp, idx_ft_all.O_info] = find_O_info(signal,labels,500,10,[],results,dscrp);
%[results, dscrp, idx_ft_all.High_Order] = find_High_Order(signal,labels,split_len,6,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_High_Order] = find_Phase_High_Order(signal,labels,split_len,3,results,dscrp);

%[results, dscrp, idx_ft_all.O_info] = find_Emergence(signal,labels,500,6,filters_PLV{1,1},results,dscrp,'_delta_band');
%[results, dscrp, idx_ft_all.O_info] = find_Emergence(signal,labels,500,6,filters_PLV{1,2},results,dscrp,'_theta_band');
%[results, dscrp, idx_ft_all.O_info] = find_Emergence(signal,labels,500,6,filters_PLV{1,3},results,dscrp,'_alpha_band');
%[results, dscrp, idx_ft_all.O_info] = find_Emergence(signal,labels,500,6,filters_PLV{1,4},results,dscrp,'_beta_band');

%[results, dscrp, idx_ft_all.Phase_Locking_Value] = find_PLV(signal,labels,256,filters_PLV,results,dscrp);
%[results, dscrp, idx_ft_all.Frequency_Mixing] = find_frequency_mixing_single_ch(signal,labels,split_len,results,dscrp,cfg);
%[results, dscrp, idx_ft_all.Frequency_Mixing] = find_frequency_mixing_cross_ch_quad(signal,labels,split_len,results,dscrp,cfg);

% single channel features
% [results, dscrp, idx_ft_all.Single_channel] = find_JL_ft(signal,f,labels,results,dscrp,filters_PLV);
% [results, dscrp, idx_ft_all.Spectral_Slope] = find_Spectral_Slope(signal,f,labels,[0.1 30],'/usr/local/bin/python3',1,results,dscrp);
% [results, dscrp, idx_ft_all.Complexity] = find_Complexity(signal,labels,results,dscrp);
% average channel
[results, dscrp, idx_ft_all.JL_features] = find_JL_ft_avg(signal,f,labels,results,dscrp,filters_PLV);
%[results, dscrp, idx_ft_all.Spectral_Slope] = find_Spectral_Slope_avg(signal,f,labels,[0.1 30],'/usr/local/bin/python3',1,results,dscrp);
%[results, dscrp, idx_ft_all.Complexity] = find_Complexity_avg(signal,labels,results,dscrp);


dscrp = transpose(dscrp);