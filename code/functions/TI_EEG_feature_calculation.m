function [results, dscrp, idx_ft_all] = TI_EEG_feature_calculation(y,ch_labels,~,f,filters_1,filters_2,cfg)
% y is the signal
% ch_label is the label for channels
% filters_1 and filters_2 are optional input filters
% cfg is configuration file for frequency mixing

split_len = 1000; % we could take the signal into even smaller parts and compute the average
signal = y;
labels = ch_labels;
results = [];
dscrp = {};
idx_ft_all = [];

%% Feature calculation
% all function follows the format of find_XX(signal, labels, split_len(you could customize), ...., results, dscrp,...);
% some might need some extra input such as filters... see more on each function

[results, dscrp, idx_ft_all.Pearson] = find_Pearson(signal,labels,split_len,results,dscrp);
[results, dscrp, idx_ft_all.Spearman] = find_Spearman(signal,labels,128,results,dscrp);
%[results, dscrp, idx_ft_all.Distcorr] = find_Distcorr(signal,labels,split_len,results,dscrp); % SW 1st one uncommented... runs
%[results, dscrp, idx_ft_all.Bandpass_Spearman] = find_bp_Spearman(signal,labels,256*6,results,dscrp,filters_2); % SW 2nd
%[results, dscrp, idx_ft_all.Bandpass_Distcorr] = find_bp_Distcorr(signal,labels,256*6,results,dscrp,filters_2);
%[results, dscrp, idx_ft_all.Cross_Correlation] = find_XCOR(signal,labels,256,results,dscrp);
%[results, dscrp, idx_ft_all.Cross_Covariance] = find_XCOV(signal,labels,256,results,dscrp);
[results, dscrp, idx_ft_all.PSDs] = find_PSDs(signal,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.FFT] = find_FFT(signal,labels,256,results,dscrp);
[results, dscrp, idx_ft_all.Power_Band] = find_PMTM(signal,f,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.Power_Band] = find_PMTM_low_high(signal,labels,split_len,results,dscrp);
[results, dscrp, idx_ft_all.Coherence] = find_CXY(signal,f,labels,3*split_len,results,dscrp);
%[results, dscrp, idx_ft_all.DTF_PDC] = find_DTF(signal,f,labels,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_Amplitude_Coupling] = find_PAC(signal,filters_2,labels,512,results,dscrp);
% [results, dscrp, idx_ft_all.Synchronization_Likelihood] = find_SL(signal,labels,split_len,results,dscrp);
% [results, dscrp, idx_ft_all.Bispectrum] = find_Bispectrum(signal,f_s,labels,split_len,results,dscrp);
% [results, dscrp, idx_ft_all.Cross_Bispectrum] = find_Cross_Bispectrum(signal,f_s,labels,split_len,results,dscrp);
% [results, dscrp, idx_ft_all.Transfer_Entropy] = find_TE_final(signal,labels,split_len,results,dscrp);
% [results, dscrp, idx_ft_all.Transfer_Entropy_pink_reduced] = find_TE_pink_reduced(signal,labels,256,results,dscrp);
% [results, dscrp, idx_ft_all.Phase_Transfer_Entropy] = find_phase_TE(signal,labels,256*6,results,dscrp);
% [results, dscrp, idx_ft_all.Phase_Transfer_Entropy] = find_BP_phase_TE(signal,labels,256*6,results,dscrp,filters_2);
% [results, dscrp, idx_ft_all.Phase_Transfer_Entropy] = find_BP_phase_TE_within_ch(signal,labels,256*6,results,dscrp,filters_2);
% [results, dscrp, idx_ft_all.Wavelet_Transfer_Entropy] = find_wavelet_TE(signal,labels,256,results,dscrp);
% [results, dscrp, idx_ft_all.Wavelet_Transfer_Entropy] = find_wavelet_TE_within_ch(signal,labels,256,results,dscrp);
% [results, dscrp, idx_ft_all.BP_Transfer_Entropy] = find_TE_BP(signal,labels,256,results,dscrp,filters_1);
% [results, dscrp, idx_ft_all.BP_Transfer_Entropy] = find_BP_TE_within_ch(signal,labels,256*6,results,dscrp,filters_2);  % SW changed last argument from filters_1 to filters_2, because filters_1 isn't actually defined in the demo_code
% [results, dscrp, idx_ft_all.Mutual_Information] = find_MI(signal,labels,split_len,results,dscrp);
% [results, dscrp, idx_ft_all.Temporal_Mutual_Information] = find_TMI(signal,labels,128,results,dscrp);
% [results, dscrp, idx_ft_all.Auto_Mutual_Information] = find_AMI(signal,labels,split_len/2,results,dscrp);
% [results, dscrp, idx_ft_all.Spectral_Entropy] = find_SE(signal,f_s,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.O_info] = find_O_info(signal,labels,6*split_len,6,filters_1,results,dscrp);
%[results, dscrp, idx_ft_all.High_Order] = find_High_Order(signal,labels,split_len,4,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_High_Order] = find_Phase_High_Order(signal,labels,split_len,3,results,dscrp);
%[results, dscrp, idx_ft_all.Phase_Locking_Value] = find_PLV(signal,labels,split_len,filters_2,results,dscrp);
%[results, dscrp, idx_ft_all.Modulation_Index] = find_Mod_index(signal,labels,split_len,results,dscrp,filters_2);
% [results, dscrp, idx_ft_all.Single_channel] = find_JL_ft(signal,f,labels,results,dscrp,filters_2); % single channel features
% [results, dscrp, idx_ft_all.Spectral_Slope] = find_Spectral_Slope_avg(signal,f,labels,[0.1 30],'/usr/local/bin/python3',1,results,dscrp);
% [results, dscrp, idx_ft_all.Complexity] = find_Complexity_avg(signal,labels,results,dscrp);

dscrp = transpose(dscrp);