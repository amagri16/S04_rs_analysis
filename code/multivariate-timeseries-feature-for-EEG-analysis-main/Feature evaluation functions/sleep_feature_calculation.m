function [results, dscrp, idx_ft_all] = sleep_feature_calculation(y,ch_labels,discard,f,filters,filters_PLV,cfg)
f_s = f;
split_len = 256*6;
filters = filters;
overlap = 0.5;
signal = y;
labels = ch_labels;
results = [];
dscrp = {};
idx_ft_all = [];

%% Feature calculation
%[results, dscrp, idx_ft_all.Pearson] = find_Pearson(signal,labels,128,results,dscrp);
%[results, dscrp, idx_ft_all.Spearman] = find_Spearman(signal,labels,128,results,dscrp);
%[results, dscrp, idx_ft_all.Distcorr] = find_Distcorr(signal,labels,128,results,dscrp);
%[results, dscrp, idx_ft_all.bp_Spearman] = find_bp_Spearman(signal,labels,128,filters,results,dscrp);
%[results, dscrp, idx_ft_all.bp_Distcorr] = find_bp_Distcorr(signal,labels,128,filters,results,dscrp);
%[results, dscrp, idx_ft_all.Cross_Correlation] = find_XCOR(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.Cross_Covariance] = find_XCOV(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.PSDs] = find_PSDs(signal,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.FFT] = find_FFT(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.Power_Band] = find_PMTM(signal,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.Power_Band] = find_PMTM_low_high(signal,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.Coherence] = find_CXY(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.DTF_PDC] = find_DTF(signal,f_s,labels,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_Amplitude_Coupling] = find_PAC(signal,filters_PLV,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.Synchronization_Likelihood] = find_SL(signal,labels,256*6,results,dscrp);
%[results, dscrp, idx_ft_all.Bispectrum] = find_Bispectrum(signal,f_s,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.Cross_Bispectrum] = find_Cross_Bispectrum(signal,f_s,labels,split_len,results,dscrp);

[results, dscrp, idx_ft_all.Transfer_Entropy] = find_TE_final(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.Transfer_Entropy_pink_reduced] = find_TE_pink_reduced(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_Transfer_Entropy] = find_phase_TE(signal,labels,256*6,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_Transfer_Entropy] = find_BP_phase_TE(signal,labels,256*6,results,dscrp,filters_PLV);
%[results, dscrp, idx_ft_all.Phase_Transfer_Entropy] = find_BP_phase_TE_within_ch(signal,labels,256*6,results,dscrp,filters_PLV);
% 
%[results, dscrp, idx_ft_all.Wavelet_Transfer_Entropy] = find_wavelet_TE(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.Wavelet_Transfer_Entropy] = find_wavelet_TE_within_ch(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.BP_Transfer_Entropy] = find_TE_BP(signal,labels,256,results,dscrp,filters);
%[results, dscrp, idx_ft_all.BP_Transfer_Entropy] = find_BP_TE_within_ch(signal,labels,256*6,results,dscrp,filters);

%[results, dscrp, idx_ft_all.Mutual_Information] = find_MI(signal,labels,128,results,dscrp);
%[results, dscrp, idx_ft_all.Temporal_Mutual_Information] = find_TMI(signal,labels,128,results,dscrp);
%[results, dscrp, idx_ft_all.Auto_Mutual_Information] = find_AMI(signal,labels,128,results,dscrp);
 
%[results, dscrp, idx_ft_all.SE] = find_SE(signal,f_s,labels,split_len,results,dscrp);

%[results, dscrp, idx_ft_all.O_info] = find_O_info(signal,labels,split_len,3,results,dscrp);
%[results, dscrp, idx_ft_all.High_Order] = find_High_Order(signal,labels,split_len,3,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_High_Order] = find_Phase_High_Order(signal,labels,split_len,3,results,dscrp);

%[results, dscrp, idx_ft_all.Phase_Locking_Value] = find_PLV(signal,labels,256,filters_PLV,results,dscrp);

%[results, dscrp, idx_ft_all.Frequency_Mixing] = find_frequency_mixing_single_ch(signal,labels,split_len,results,dscrp,cfg);
%[results, dscrp, idx_ft_all.Frequency_Mixing] = find_frequency_mixing_cross_ch_quad(signal,labels,split_len,results,dscrp,cfg);
dscrp = transpose(dscrp);