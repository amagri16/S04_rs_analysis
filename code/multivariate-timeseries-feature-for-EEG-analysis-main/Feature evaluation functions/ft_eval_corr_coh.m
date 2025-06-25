function sampleObj_coll = ft_eval_corr_coh(sampleObj_coll,filters,filters_PLV,cfg)

% This function is used to evaluate the features and turn it into a feature
% vector with corresponding descriptions
%
% Author: Junheng Li/Tianyu Wei

%% Log of code

% 29/01/2023, code started base on ft_eval.m


% to do list:
    % use xcorr to find the max value and its lag, average of all the
    % values
%% Feature evaluations

num_samps = length(sampleObj_coll);        % Total number of samples in the collection

% Enables parallel computation 
parfor i = 1:num_samps %iterating through the number of samples 
    samp_now = sampleObj_coll{i}; %extract the ith sample structure

    y = samp_now.data_EEGs;
    ch_labels = samp_now.Channel_labels;
    Fs = samp_now.Fs;
    discard_info = samp_now.discard_info;
    [ft,ft_dscrp,idx] = MCI_feature_calculation(y,ch_labels,discard_info,Fs,filters,filters_PLV,cfg);        % Core feature evaluation
    samp_now.ft_vec = ft;
    samp_now.ft_idx = idx;
%     if sum(isnan(ft))>0
%         samp_now.ft_quality = 0;
%         disp(['Sample No.',num2str(i),' from Subject No.',num2str(samp_now.Sbj_id),' have NaN feature value'])
%         continue
%     end
    samp_now.ft_quality = 1;
    samp_now.ft_dscrp = ft_dscrp;
    %If all is okay,siognifiy it by the ft_quality variable = 1 and
    %print out that ft-eval went well 
    disp(['Sample No.',num2str(i),' from Subject No.',num2str(samp_now.Sbj_id),' has been successfully calculated'])
    sampleObj_coll{i} = samp_now; %input the extracted features and corresponding iformation back into teh samples form the subject structure 
    end
end
