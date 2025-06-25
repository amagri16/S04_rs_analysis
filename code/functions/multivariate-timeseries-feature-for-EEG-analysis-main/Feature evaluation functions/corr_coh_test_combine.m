function [results, dscrp] = corr_coh_test_combine(y,ch_labels,discard,f)
f_s = f;
split_len = 256*4;
nfft = 0:0.05:30;
delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [13,30];
signal = y;
labels = ch_labels;
discard_info = discard;
ava_ch_num = size(signal);
ava_ch_index = 1;
results = [];
dscrp = {};

for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};

        % assign the description of the features
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_pearson';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_spearman';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_half_max_xorr_lag_width';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_half_max_xcov_lag_width';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_coherence';

        %         if discard_info(m) | discard_info(n)
        %             results = [results nan nan nan nan nan nan nan nan];
        %         else


        % creates container for feature value from each epoch
        pearson_all = [];
        spearman_all = [];
        spearman_psd_all = [];
        delta_spearman_psd_all = [];
        theta_spearman_psd_all = [];
        alpha_spearman_psd_all = [];
        beta_spearman_psd_all = [];

        mean_delta_csd_all = [];
        mean_theta_csd_all = [];
        mean_alpha_csd_all = [];
        mean_beta_csd_all = [];

        half_max_xorr_lag_width_all = [];
        half_max_xcov_lag_width_all = [];

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
            % add feature analysis for each signal
             % pearson
            sig1_now = sig1_split{k};
            sig2_now = sig2_split{k};
            sig_size = split_len; 

            R_pearson_global = corrcoef(sig1_now,sig2_now);
            pearson_all = [pearson_all R_pearson_global(2,1)];

            spearman_corr = corr(sig1_now, sig2_now, 'Type', 'Spearman');
            spearman_all = [spearman_all spearman_corr];

            % global xcorr
            [R,lags] = xcorr(sig1_now,sig2_now);
            [corr_max,corr_max_index] = max(R);
            max_lag = lags(corr_max_index);
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
            half_max = 0.5*corr_max;
            % all the index that is larger than half of the max
            corr_half_max_lags = (R>half_max);
            [labeledX, ~] = bwlabel(corr_half_max_lags);
            props = regionprops(labeledX, 'Area', 'PixelList');
            regionLengths = [props.Area];
            half_lag = find(labeledX==labeledX(corr_max_index));
            half_lag_len = size(half_lag,1);
            half_max_xcov_lag_width_all = [half_max_xcov_lag_width_all half_lag_len];


            % Compute the power spectral density of the signal using Welch's
            % method
            [psd1,~] = pwelch(sig1_now, [], [], 2560, 256);
            [psd2,~] = pwelch(sig2_now, [], [], 2560, 256);

            % compute pearson
            spearman_psd = corr(psd1,psd2, 'Type', 'Spearman');
            spearman_psd_all = [spearman_psd_all spearman_psd];
            
            delta_spearman_psd = corr(psd1([11:1:41]),psd2([11:1:41]), 'Type', 'Spearman');
            theta_spearman_psd = corr(psd1([41:1:81]),psd2([41:1:81]), 'Type', 'Spearman');
            alpha_spearman_psd = corr(psd1([81:1:121]),psd2([81:1:121]), 'Type', 'Spearman');
            beta_spearman_psd = corr(psd1([131:1:301]),psd2([131:1:301]), 'Type', 'Spearman');

            delta_spearman_psd_all = [delta_spearman_psd_all delta_spearman_psd];
            theta_spearman_psd_all = [theta_spearman_psd_all theta_spearman_psd];
            alpha_spearman_psd_all = [alpha_spearman_psd_all alpha_spearman_psd];
            beta_spearman_psd_all = [beta_spearman_psd_all beta_spearman_psd];


            % compute coherence

            [csd, freqs] = mscohere(sig1_now, sig2_now,[],[], nfft, f_s);
            delta_idx = freqs >= delta_bd(1) & freqs <= delta_bd(2);
            theta_idx = freqs >= theta_bd(1) & freqs <= theta_bd(2);
            alpha_idx = freqs >= alpha_bd(1) & freqs <= alpha_bd(2);
            beta_idx = freqs >= beta_bd(1) & freqs <= beta_bd(2);
            csd_delta = sum(csd(delta_idx));
            csd_theta = sum(csd(theta_idx));
            csd_alpha = sum(csd(alpha_idx));
            csd_beta = sum(csd(beta_idx));
            
            mean_delta_csd_all = [mean_delta_csd_all csd_delta];
            mean_theta_csd_all = [mean_theta_csd_all csd_theta];
            mean_alpha_csd_all = [mean_alpha_csd_all csd_alpha];
            mean_beta_csd_all = [mean_beta_csd_all csd_beta];
           
        end
        results = [results mean(pearson_all)];
        results = [results mean(spearman_all)];
        results = [results mean(half_max_xorr_lag_width_all)];
        results = [results mean(half_max_xcov_lag_width_all)];
        results = [results mean(spearman_psd_all)];
        results = [results mean(delta_spearman_psd_all)];
        results = [results mean(theta_spearman_psd_all)];
        results = [results mean(alpha_spearman_psd_all)];
        results = [results mean(beta_spearman_psd_all)];

        results = [results mean(mean_delta_csd_all)];
        results = [results mean(mean_theta_csd_all)];
        results = [results mean(mean_alpha_csd_all)];
        results = [results mean(mean_beta_csd_all)];
       
        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
dscrp = transpose(dscrp);
