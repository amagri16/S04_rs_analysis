function [results, dscrp] = corr_coh_test_combine_ensembled(y,ch_labels,discard,f)
f_s = f;
split_len = 256*4;
nfft = 0:0.05:30;
delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [12,30];
beta_bd_1 = [12,18];
beta_bd_2 = [18,24];
beta_bd_3 = [24,30];
%
% make change to this
% psd_y_norm = psd_y ./ sum(psd_y);    % Normalize to total power
%
signal = y;
labels = ch_labels;
discard_info = discard;
ava_ch_num = size(signal);
ava_ch_index = 1;
results = [];
dscrp = {};
d_delta = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',delta_bd(1),'CutoffFrequency2',delta_bd(2),'SampleRate',f_s);
d_theta = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',theta_bd(1),'CutoffFrequency2',theta_bd(2),'SampleRate',f_s);
d_alpha = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',alpha_bd(1),'CutoffFrequency2',alpha_bd(2),'SampleRate',f_s);
d_beta = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',beta_bd_1(1),'CutoffFrequency2',beta_bd_1(2),'SampleRate',f_s);
filters = {d_delta,d_theta,d_alpha,d_beta};


for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};

        % assign the description of the features
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_pearson';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_Spearman';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_half_max_xorr_lag_width';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_max_xorr_lag';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_half_max_xcov_lag_width';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_max_xcov_lag';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_energy';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_spectral_flatness';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_Spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_Spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_Spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_Spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_psd_1';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_psd_2';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_psd_3';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_psd_1';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_psd_2';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_psd_3';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_coherence_1';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_coherence_2';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_coherence_3';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_coherence_flatness_all_bands';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_SP_Summaries_welch_rect_area_5_1';

        pearson_all = [];
        Spearman_all = [];
        Spearman_psd_all = [];

        delta_psd1_all = [];
        theta_psd1_all = [];
        alpha_psd1_all = [];
        beta_psd1_all = [];
        beta1_psd1_all = [];
        beta2_psd1_all = [];
        beta3_psd1_all = [];

        delta_psd2_all = [];
        theta_psd2_all = [];
        alpha_psd2_all = [];
        beta_psd2_all = [];
        beta1_psd2_all = [];
        beta2_psd2_all = [];
        beta3_psd2_all = [];

        delta_Spearman_psd_all = [];
        theta_Spearman_psd_all = [];
        alpha_Spearman_psd_all = [];
        beta_Spearman_psd_all = [];
        beta1_Spearman_psd_all = [];
        beta2_Spearman_psd_all = [];
        beta3_Spearman_psd_all = [];

        mean_delta_csd_all = [];
        mean_theta_csd_all = [];
        mean_alpha_csd_all = [];
        mean_beta_csd_all = [];
        mean_beta1_csd_all = [];
        mean_beta2_csd_all = [];
        mean_beta3_csd_all = [];

        half_max_xorr_lag_width_all = [];
        max_xorr_lag_all = [];
        half_max_xcov_lag_width_all = [];
        max_xcov_lag_all = [];

        energy_sig1 = [];
        spectral_flatness_sig1 = [];
        energy_sig2 = [];
        spectral_flatness_sig2 = [];
        
        coherence_flatness = [];

        SP_Summaries_welch_rect_area_5_1_sig1 = [];
        SP_Summaries_welch_rect_area_5_1_sig2 = [];
        
        sig1_split = {};
        sig2_split = {};
        num_sig1 = floor(size(sig1,1)/split_len);
        num_sig2 = floor(size(sig2,1)/split_len);
        % split the data into sub parts
        for i=0:1:num_sig1-1
            sig1_split{i+1,1}=sig1(i*split_len+1:1:i*split_len+(split_len));
        end

        for i=0:1:num_sig2-1
            sig2_split{i+1,1}=sig2(i*split_len+1:1:i*split_len+(split_len));
        end

        sig1_split_size = size(sig1_split,1);
        sig2_split_size = size(sig2_split,1);
        min_sig_size = min(sig1_split_size,sig2_split_size);
        % traverse through the splited signals

        for k = 1:1:min_sig_size
            % add feature analysis for each signal pearson
            sig1_now = sig1_split{k};
            sig2_now = sig2_split{k};
            sig_size = split_len; 

            R_pearson_global = corrcoef(sig1_now,sig2_now);
            pearson_all = [pearson_all R_pearson_global(2,1)];

            Spearman_corr = corr(sig1_now, sig2_now, 'Type', 'Spearman');
            Spearman_all = [Spearman_all Spearman_corr];

            % global xcorr
            [R,lags] = xcorr(sig1_now,sig2_now);
            [corr_max,corr_max_index] = max(R);
            max_lag = lags(corr_max_index);
            max_xorr_lag_all = [max_xorr_lag_all; max_lag];
            half_max = 0.5*corr_max;
            % all the index that is larger than half of the max
            corr_half_max_lags = (R>half_max);
            [labeledX, ~] = bwlabel(corr_half_max_lags);
            half_lag = find(labeledX==labeledX(corr_max_index));
            half_lag_len = size(half_lag,1);

            half_max_xorr_lag_width_all = [half_max_xorr_lag_width_all half_lag_len];


            %  global xcov
            [R,lags] = xcov(sig1_now,sig2_now);
            [corr_max,corr_max_index] = max(R);
            max_lag = lags(corr_max_index);
            max_xcov_lag_all = [max_xcov_lag_all;max_lag];
            half_max = 0.5*corr_max;
            % all the index that is larger than half of the max
            corr_half_max_lags = (R>half_max);
            [labeledX, ~] = bwlabel(corr_half_max_lags);
            props = regionprops(labeledX, 'Area', 'PixelList');
            regionLengths = [props.Area];
            half_lag = find(labeledX==labeledX(corr_max_index));
            half_lag_len = size(half_lag,1);
            half_max_xcov_lag_width_all = [half_max_xcov_lag_width_all half_lag_len];


            % Compute the energy
            energy_val_sig1 = sum(sig1_now.^2);
            energy_val_sig2 = sum(sig2_now.^2);

            energy_sig1 = [energy_sig1; energy_val_sig1];
            energy_sig2 = [energy_sig2; energy_val_sig2];


            % Compute the power spectral density of the signal using Welch's
            % method

            sig1_diff = sig1_now(2:end) - sig1_now(1:end-1);
            sig2_diff = sig2_now(2:end) - sig2_now(1:end-1);
            [psd1,~] = pwelch(sig1_diff, [], [], 2560, 256);
            [psd2,~] = pwelch(sig2_diff, [], [], 2560, 256);
            
            psd1 = psd1/sum(psd1);
            psd2 = psd2/sum(psd2);

            % Compute the spectral flatness
            spectral_flatness_val_sig1 = geomean(psd1) / mean(psd1);
            spectral_flatness_val_sig2 = geomean(psd2) / mean(psd2);
            
            spectral_flatness_sig1 = [spectral_flatness_sig1; spectral_flatness_val_sig1];
            spectral_flatness_sig2 = [spectral_flatness_sig2; spectral_flatness_val_sig2];

            % compute pearson
            Spearman_psd = corr(psd1,psd2, 'Type', 'Spearman');
            Spearman_psd_all = [Spearman_psd_all Spearman_psd];
            
            delta_Spearman_psd = corr(psd1([11:1:41]),psd2([11:1:41]), 'Type', 'Spearman');
            theta_Spearman_psd = corr(psd1([41:1:81]),psd2([41:1:81]), 'Type', 'Spearman');
            alpha_Spearman_psd = corr(psd1([81:1:121]),psd2([81:1:121]), 'Type', 'Spearman');
            beta_Spearman_psd = corr(psd1([121:1:301]),psd2([121:1:301]), 'Type', 'Spearman');
            beta1_Spearman_psd = corr(psd1([121:1:181]),psd2([121:1:181]), 'Type', 'Spearman');
            beta2_Spearman_psd = corr(psd1([181:1:241]),psd2([181:1:241]), 'Type', 'Spearman');
            beta3_Spearman_psd = corr(psd1([241:1:301]),psd2([241:1:301]), 'Type', 'Spearman');

            delta_Spearman_psd_all = [delta_Spearman_psd_all delta_Spearman_psd];
            theta_Spearman_psd_all = [theta_Spearman_psd_all theta_Spearman_psd];
            alpha_Spearman_psd_all = [alpha_Spearman_psd_all alpha_Spearman_psd];
            beta_Spearman_psd_all = [beta_Spearman_psd_all beta_Spearman_psd];
            beta1_Spearman_psd_all = [beta1_Spearman_psd_all beta1_Spearman_psd];
            beta2_Spearman_psd_all = [beta2_Spearman_psd_all beta2_Spearman_psd];
            beta3_Spearman_psd_all = [beta3_Spearman_psd_all beta3_Spearman_psd]; 
            
            delta_psd1_all = [delta_psd1_all; mean(psd1([11:1:41]))];
            theta_psd1_all = [theta_psd1_all; mean(psd1([41:1:81]))];
            alpha_psd1_all = [alpha_psd1_all; mean(psd1([81:1:121]))];
            beta_psd1_all = [beta_psd1_all; mean(psd1([121:1:301]))];
            beta1_psd1_all = [beta1_psd1_all; mean(psd1([121:1:181]))];
            beta2_psd1_all = [beta2_psd1_all; mean(psd1([181:1:241]))];
            beta3_psd1_all = [beta3_psd1_all; mean(psd1([241:1:301]))];

            delta_psd2_all = [delta_psd2_all; mean(psd2([11:1:41]))];
            theta_psd2_all = [theta_psd2_all; mean(psd2([41:1:81]))];
            alpha_psd2_all = [alpha_psd2_all; mean(psd2([81:1:121]))];
            beta_psd2_all = [beta_psd2_all; mean(psd2([121:1:301]))];
            beta1_psd2_all = [beta1_psd2_all; mean(psd2([121:1:181]))];
            beta2_psd2_all = [beta2_psd2_all; mean(psd2([181:1:241]))];
            beta3_psd2_all = [beta3_psd2_all; mean(psd2([241:1:301]))];
            % compute coherence density

            [csd, freqs] = mscohere(sig1_now, sig2_now,[],[], nfft, f_s);
            delta_idx = freqs >= delta_bd(1) & freqs <= delta_bd(2);
            theta_idx = freqs >= theta_bd(1) & freqs <= theta_bd(2);
            alpha_idx = freqs >= alpha_bd(1) & freqs <= alpha_bd(2);
            beta_idx = freqs >= beta_bd(1) & freqs <= beta_bd(2);
            beta1_idx = freqs >= beta_bd_1(1) & freqs <= beta_bd_1(2);
            beta2_idx = freqs >= beta_bd_2(1) & freqs <= beta_bd_2(2);
            beta3_idx = freqs >= beta_bd_3(1) & freqs <= beta_bd_3(2);

            csd_delta = sum(csd(delta_idx));
            csd_theta = sum(csd(theta_idx));
            csd_alpha = sum(csd(alpha_idx));
            csd_beta = sum(csd(beta_idx));
            csd_beta1 = sum(csd(beta1_idx));            
            csd_beta2 = sum(csd(beta2_idx));            
            csd_beta3 = sum(csd(beta3_idx));
            
            mean_delta_csd_all = [mean_delta_csd_all csd_delta];
            mean_theta_csd_all = [mean_theta_csd_all csd_theta];
            mean_alpha_csd_all = [mean_alpha_csd_all csd_alpha];
            mean_beta_csd_all = [mean_beta_csd_all csd_beta];
            mean_beta1_csd_all = [mean_beta1_csd_all csd_beta1];
            mean_beta2_csd_all = [mean_beta2_csd_all csd_beta2];
            mean_beta3_csd_all = [mean_beta3_csd_all csd_beta3];
            coherence_flatness = [coherence_flatness geomean(csd)/mean(csd)];

            SP_Summaries_welch_rect_area_5_1_sig1 = [SP_Summaries_welch_rect_area_5_1_sig1; SP_Summaries_welch_rect_area_5_1(sig1_now)];
            SP_Summaries_welch_rect_area_5_1_sig2 = [SP_Summaries_welch_rect_area_5_1_sig2; SP_Summaries_welch_rect_area_5_1(sig2_now)];

        end
        results = [results mean(pearson_all)];
        results = [results mean(Spearman_all)];
        results = [results mean(half_max_xorr_lag_width_all)];
        results = [results mean(max_xorr_lag_all)];
        results = [results mean(half_max_xcov_lag_width_all)];
        results = [results mean(max_xcov_lag_all)];
        p_energy = corr(energy_sig1, energy_sig2,'Type', 'Spearman');
        results = [results p_energy];

        p_spectral_flatness = corr(spectral_flatness_sig1, spectral_flatness_sig2,'Type', 'Spearman');
        results = [results p_spectral_flatness];

        results = [results mean(Spearman_psd_all)];
        results = [results mean(delta_Spearman_psd_all)];
        results = [results mean(theta_Spearman_psd_all)];
        results = [results mean(alpha_Spearman_psd_all)];
        results = [results mean(beta_Spearman_psd_all)];
        results = [results mean(beta1_Spearman_psd_all)];
        results = [results mean(beta2_Spearman_psd_all)];
        results = [results mean(beta3_Spearman_psd_all)];

        delta_psd = corr(delta_psd1_all,delta_psd2_all,'Type', 'Spearman');
        theta_psd = corr(theta_psd1_all,theta_psd2_all,'Type', 'Spearman');
        alpha_psd = corr(alpha_psd1_all,alpha_psd2_all,'Type', 'Spearman');
        beta_psd = corr(beta_psd1_all,beta_psd2_all,'Type', 'Spearman');
        beta1_psd = corr(beta1_psd1_all,beta1_psd2_all,'Type', 'Spearman');
        beta2_psd = corr(beta2_psd1_all,beta2_psd2_all,'Type', 'Spearman');
        beta3_psd = corr(beta3_psd1_all,beta3_psd2_all,'Type', 'Spearman');

        p_SP_Summaries_welch_rect_area_5_1 = corr(SP_Summaries_welch_rect_area_5_1_sig1,...
            SP_Summaries_welch_rect_area_5_1_sig2,'Type', 'Spearman');

        results = [results delta_psd];
        results = [results theta_psd];
        results = [results alpha_psd];
        results = [results beta_psd];
        results = [results beta1_psd];
        results = [results beta2_psd];
        results = [results beta3_psd];

        results = [results mean(mean_delta_csd_all)];
        results = [results mean(mean_theta_csd_all)];
        results = [results mean(mean_alpha_csd_all)];
        results = [results mean(mean_beta_csd_all)];
        results = [results mean(mean_beta1_csd_all)];
        results = [results mean(mean_beta2_csd_all)];
        results = [results mean(mean_beta3_csd_all)];
        results = [results mean(coherence_flatness)];
        results = [results p_SP_Summaries_welch_rect_area_5_1]; 

        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end

[results, dscrp] = find_Pearson(signal,results,dscrp);
[results, dscrp] = find_Spearman(signal,results,dscrp);
[results, dscrp] = find_XCOR(signal,results,dscrp);
[results, dscrp] = find_XCOV(signal,results,dscrp);
[results, dscrp] = find_Energy(signal,results,dscrp);

[results, dscrp] = find_PSDs(signal,results,dscrp);
[results, dscrp] = find_CXY(signal,results,dscrp);

[results, dscrp] = find_DTF(signal,f_s,labels,results,dscrp);
[results, dscrp] = find_PAC(signal,filters,labels,split_len,results,dscrp);
[results, dscrp] = find_SL(signal,split_len,labels,results,dscrp);
dscrp = transpose(dscrp);