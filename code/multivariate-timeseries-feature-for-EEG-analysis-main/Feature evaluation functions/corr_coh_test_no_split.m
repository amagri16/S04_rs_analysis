function [results, dscrp] = corr_coh_test_no_split(y,ch_labels,discard,f)
f_s = f;
nfft = 0:0.1:30;
delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [13,30];
signal = y;
labels = ch_labels;
discard_info = discard;
ava_ch_num = size(signal);
results = [];
dscrp = {};

for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};
        sig_size = size(sig1,1);

        % assign the description of the features

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_pearson';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_half_max_xorr_lag_width';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_half_max_xcov_lag_width';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_energy';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_spectral_flatness';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_kurtosis_psd_result';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_pearson_psd';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_delta_csd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_theta_csd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_alpha_csd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_beta_csd';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CO_trev_1_num';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_SP_Summaries_welch_rect_area_5_1';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_SP_Summaries_welch_rect_centroid';



        if discard_info(m) | discard_info(n) | isempty(sig1) | isempty(sig2)
            results = [results ...
                nan nan nan nan nan nan nan...
                nan nan nan nan nan nan nan ];
        else

            % pearson
            R_pearson_global = corrcoef(sig1,sig2);
            results = [results R_pearson_global(2,1)];


            % global xcorr
            R = xcorr(sig1,sig2);
            [corr_max,corr_max_index] = max(R);
            max_lag = corr_max_index-sig_size-1;
            half_max = 0.5*corr_max;
            % all the index that is larger than half of the max
            corr_half_max_lags = (R>half_max);
            [labeledX, ~] = bwlabel(corr_half_max_lags);
            props = regionprops(labeledX, 'Area', 'PixelList');
            regionLengths = [props.Area];
            half_lag = find(labeledX==labeledX(corr_max_index));
            half_lag_len = size(half_lag,1);

            results = [results half_lag_len];


            %  global xcov
            R = xcov(sig1,sig2);
            [corr_max,corr_max_index] = max(R);
            max_lag = corr_max_index-sig_size-1;
            half_max = 0.5*corr_max;
            % all the index that is larger than half of the max
            corr_half_max_lags = (R>half_max);
            [labeledX, ~] = bwlabel(corr_half_max_lags);
            props = regionprops(labeledX, 'Area', 'PixelList');
            regionLengths = [props.Area];
            half_lag = find(labeledX==labeledX(corr_max_index));
            half_lag_len = size(half_lag,1);
            results = [results half_lag_len];


            % Compute the energy
            energy_val_sig1 = sum(sig1.^2);
            energy_val_sig2 = sum(sig2.^2);
            energy =  mean(energy_val_sig1-energy_val_sig2);


            % Compute the power spectral density of the signal using Welch's
            % method
            [psd1,~] = pwelch(sig1, [], [], 2560, 256);
            [psd2,~] = pwelch(sig2, [], [], 2560, 256);


            % Compute the spectral flatness
            spectral_flatness_val_sig1 = geomean(psd1) / mean(psd1);
            spectral_flatness_val_sig2 = geomean(psd2) / mean(psd2);
            spectral_flatness =  mean(spectral_flatness_val_sig1-spectral_flatness_val_sig2);


            % Compute the kurtosis
            kurtosis_psd_sig1 = kurtosis(psd1);
            kurtosis_psd_sig2 = kurtosis(psd2);
            kurtosis_psd_result =  mean(kurtosis_psd_sig1-kurtosis_psd_sig2);


            % compute pearson
            pearson_psd = corrcoef(psd1,psd2);
            results = [results energy spectral_flatness kurtosis_psd_result pearson_psd(2,1)];



            window_size = length(sig1); % Window size for segmenting the signals
            overlap_ratio = 0.5; % Ratio of overlap between adjacent segments

            [csd, freqs] = cpsd(sig1, sig2, window_size, window_size*overlap_ratio, nfft, f_s);
            delta_idx = freqs >= delta_bd(1) & freqs <= delta_bd(2);
            theta_idx = freqs >= theta_bd(1) & freqs <= theta_bd(2);
            alpha_idx = freqs >= alpha_bd(1) & freqs <= alpha_bd(2);
            beta_idx = freqs >= beta_bd(1) & freqs <= beta_bd(2);
            csd_delta = mean(abs(csd(delta_idx)));
            csd_theta = mean(abs(csd(theta_idx)));
            csd_alpha = mean(abs(csd(alpha_idx)));
            csd_beta = mean(abs(csd(beta_idx)));
            results = [results csd_delta csd_theta csd_alpha csd_beta];



            % catch 22 features

            [ft_values_sig1,ft_name_sig1] = catch22_all(sig1); % run catch22-master/wrap_Matlab/mexAll.m before executing this line
            [ft_values_sig2,~] = catch22_all(sig2);

            %         for i = [6,16,21]
            %             ft_name_now = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_'+ft_name_sig1{1,i};
            %             dscrp{end+1} = ft_name_now;
            %         end
            catch22_result = [];

            for i = [6,16,21]
                catch22_result = [catch22_result ft_values_sig1(i)-ft_values_sig2(i)];
            end
            results = [results catch22_result];
        end
    end
end
dscrp = transpose(dscrp);