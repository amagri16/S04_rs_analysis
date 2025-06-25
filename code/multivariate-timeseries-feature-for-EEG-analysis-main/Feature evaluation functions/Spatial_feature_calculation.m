function [results, dscrp, idx_ft_all] = Spatial_feature_calculation(y,ch_labels,f,filters_1,filters_2)
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
% some might need some extra input such as filters to analyze specific band but one edit that in each function... see more on each function

%% Bivariate analysis

%[results, dscrp, idx_ft_all.Pearson] = find_Pearson(signal,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.Spearman] = find_Spearman(signal,labels,split_len,results,dscrp); 
%[results, dscrp, idx_ft_all.Distcorr] = find_Distcorr(signal,labels,split_len,results,dscrp);
% [results, dscrp, idx_ft_all.Cross_Correlation] = find_XCOR(signal,labels,split_len,results,dscrp);

%[results, dscrp, idx_ft_all.Cross_Covariance] =
%find_XCOV(signal,labels,split_len,results,dscrp); i am commenting this
%line because find_XCOV calls bwlabel, which is part of the Image
%Processing Toolbox (IPT), I don't have the licence to install it

% [results, dscrp, idx_ft_all.FFT] = find_FFT(signal,labels,split_len,results,dscrp);
% [results, dscrp, idx_ft_all.Power_Band] = find_PMTM(signal,f,labels,split_len,results,dscrp);
%[results, dscrp, idx_ft_all.Coherence] = find_CXY(signal,f,labels,3*split_len,results,dscrp);
% [results, dscrp, idx_ft_all.Synchronization_Likelihood] = find_SL(signal,labels,split_len,results,dscrp);
% [results, dscrp, idx_ft_all.DTF_PDC] = find_DTF(signal,f,labels,results,dscrp);
% [results, dscrp, idx_ft_all.Phase_Locking_Value] = find_PLV(signal,labels,split_len,filters_2,results,dscrp);

%% Bivariate cross-frequency coupling

% [results, dscrp, idx_ft_all.Cross_Bispectrum] = find_Cross_Bispectrum(signal,f,labels,split_len,results,dscrp);
% [results, dscrp, idx_ft_all.Modulation_Index] = find_Mod_index(signal,labels,split_len,results,dscrp,filters_2);

%% Bivariate information theory based metrics

% [results, dscrp, idx_ft_all.Transfer_Entropy] = find_TE_final(signal,labels,split_len,results,dscrp); % compute transfer entropy
% [results, dscrp, idx_ft_all.Mutual_Information] = find_MI(signal,labels,split_len,results,dscrp); % mutual information
% [results, dscrp, idx_ft_all.Temporal_Mutual_Information] = find_TMI(signal,labels,split_len,results,dscrp); % temporally resolved MI

%% Higher-order interactions
% two different packages for O-information computation
% [results, dscrp, idx_ft_all.O_info] = find_O_info(signal,labels,split_len,6,[],results,dscrp,''); 
% [results, dscrp, idx_ft_all.High_Order] = find_High_Order(signal,labels,split_len,4,results,dscrp);

%% Single channel analysis

[results, dscrp, idx_ft_all.Single_channel] = find_JL_ft(signal,f,labels,results,dscrp,filters_1); % single channel features from Junheng's work
% [results, dscrp, idx_ft_all.Added_fts] = JL_paper_version_added_fts(signal,f,labels,[0.1 30],'/usr/local/bin/python3.10',1,results,dscrp); % single channel features from Junheng's work

% [results, dscrp, idx_ft_all.Spectral_Slope] = find_Spectral_Slope_avg(signal,f,labels,[0.1 30],'/usr/local/bin/python3',1,results,dscrp);
% [results, dscrp, idx_ft_all.Complexity] = find_Complexity_avg(signal,labels,results,dscrp);
% [results, dscrp, idx_ft_all.Bispectrum] = find_Bispectrum(signal,f,labels,split_len,results,dscrp);

dscrp = transpose(dscrp);
